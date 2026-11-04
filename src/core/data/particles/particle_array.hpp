#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP


#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/utilities/range/range.hpp"

#include "particle.hpp"

namespace PHARE::core
{
// scratch of the domain bag sort, owned by the caller and reused across calls
template<typename Particle_t>
struct ParticleSortScratch
{
    std::vector<Particle_t> particles; // full layout: scatter target (swapped with the array)
    std::vector<int> keys, counts, csr, ends, bufKeys;
    std::vector<Particle_t> buf;
    std::vector<Particle_t> leavers; // particles whose key is the tail bucket
    std::vector<int> vacated;        // slots vacated by the movers, ascending
};

/** Bag layout of a sorted array: capacity of a row at a full layout = count + bagSlots for
 * non-empty rows, 0 for empty ones; a full row borrows a slot from the nearest row with a free
 * slot within bagRadius rows (cascade), else the whole array is laid out again.
 * The slack is small on purpose: gaps of more than a couple of particles between rows break the
 * hardware prefetch stream of the kernels that walk the rows (measured: 5 % slack costs ~10 %
 * of the 2D interp-1 push, 20 % slack ~50 %).
 */
inline constexpr int bagSlots  = 2;
inline constexpr int bagRadius = 256;

/** Particle array with a cell map, and an optional sorted mode.
 *
 * Sorted mode: csr_ (size nCells+2) over box_, key(cell) = row-major local flat index, last
 * index fastest (field memory order); key nCells is the tail bucket (leavers, see sortBags).
 * Rows are bags (sortBags: row k = [csr_[k], ends_[k]), holes up to csr_[k + 1]) or, after a
 * compaction, a packed CSR (row k = [csr_[k], csr_[k + 1])).
 * csrValid_: particles [0, csrPrefix_) sit in their key row; particles appended after the sort
 * ([csrPrefix_, size)) are an unsorted tail that box queries scan linearly. Appends keep the
 * CSR, any other mutation through this API drops it. In-place iCell writes through operator[]
 * are NOT tracked: the caller that pushes in place calls invalidateCSR() first (ion_updater).
 * mapValid_: the cell map indexes every particle. dropMap() stops maintaining it (sorted
 * mode); a query that needs it while csr is invalid rebuilds it lazily. An array that is never
 * sorted keeps mapValid_ true and csrValid_ false, and behaves as a plain cell-mapped array.
 *
 * Compaction contract of a bag layout: size() and operator[] are the raw storage, holes
 * included. Code that walks [0, size()) either walks the live runs (forEachLiveRun, the rows of
 * csr()/rowEnds()) or calls compact() first; begin()/end(), vector(), resize(), operator== and
 * the cell-map queries compact on their own.
 */
template<std::size_t dim>
class ParticleArray
{
public:
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = dim;
    using This                          = ParticleArray<dim>;
    using Particle_t                    = Particle<dim>;
    using Vector                        = std::vector<Particle_t>;
    using SortScratch                   = ParticleSortScratch<Particle_t>;

private:
    using CellMap_t   = CellMap<dim, int>;
    using IndexRange_ = IndexRange<This>;


public:
    using value_type     = Particle_t;
    using box_t          = Box<int, dim>;
    using iterator       = Vector::iterator;
    using const_iterator = Vector::const_iterator;



public:
    ParticleArray(box_t box)
        : box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
        initShape_();
    }

    ParticleArray(box_t box, std::size_t size)
        : particles_(size)
        , box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
        initShape_();
    }

    ParticleArray(ParticleArray const& from)            = default;
    ParticleArray(ParticleArray&& from)                 = default;
    ParticleArray& operator=(ParticleArray&& from)      = default;
    ParticleArray& operator=(ParticleArray const& from) = default;

    NO_DISCARD std::size_t size() const { return particles_.size(); }
    NO_DISCARD std::size_t capacity() const { return particles_.capacity(); }

    void clear()
    {
        particles_.clear();
        cellMap_.clear();
        mapValid_ = true;
        csrValid_ = false;
        bagged_   = false;
    }
    void reserve(std::size_t newSize) { return particles_.reserve(newSize); }
    void resize(std::size_t newSize)
    {
        compact();
        csrValid_ = false;
        return particles_.resize(newSize);
    }

    NO_DISCARD auto const& operator[](std::size_t i) const { return particles_[i]; }
    NO_DISCARD auto& operator[](std::size_t i) { return particles_[i]; }

    NO_DISCARD bool operator==(ParticleArray<dim> const& that) const
    {
        compact();
        that.compact();
        return (this->particles_ == that.particles_);
    }

    // bags: iterators see the dense array, the holes are compacted away
    // first. size() and operator[] are the raw storage.
    NO_DISCARD auto begin() const { return compact(), particles_.begin(); }
    NO_DISCARD auto begin() { return compact(), particles_.begin(); }

    NO_DISCARD auto end() const { return compact(), particles_.end(); }
    NO_DISCARD auto end() { return compact(), particles_.end(); }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        assert(!bagged_); // position comes from begin()
        csrValid_ = false;
        particles_.insert(position, first, last);
    }

    NO_DISCARD auto back() { return particles_.back(); }
    NO_DISCARD auto front() { return particles_.front(); }


    auto erase(IndexRange_ range)
    {
        assert(!bagged_); // the range comes from a dense array (partition)
        if (mapValid_)
            cellMap_.erase(range); // erases the items through erase(iterator, iterator)
        else
            erase(range.begin(), range.end());
    }

    iterator erase(iterator first, iterator last)
    {
        // should we erase particles indexes associated with these iterators from the cellmap?
        // probably it does not matter if not. The reason is that
        // particles erased from the particlearray are so because they left
        // the patch cells to an outside cell.
        // But in principle that cell will never be accessed because it is outside the patch.
        // The only thing "bad" if these indexes are not deleted is that the
        // size of the cellmap becomes unequal to the size of the particleArray.
        // but  ¯\_(ツ)_/¯
        if (static_cast<std::size_t>(first - particles_.begin()) < csrPrefix_)
            csrValid_ = false;
        return particles_.erase(first, last);
    }


    Particle_t& emplace_back()
    {
        auto& part = particles_.emplace_back();
        if (mapValid_)
            cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }



    Particle_t& emplace_back(Particle_t&& p)
    {
        auto& part = particles_.emplace_back(std::forward<Particle_t>(p));
        if (mapValid_)
            cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    void push_back(Particle_t const& p)
    {
        particles_.push_back(p);
        if (mapValid_)
            cellMap_.add(particles_, particles_.size() - 1);
    }

    void push_back(Particle_t&& p)
    {
        particles_.push_back(std::forward<Particle_t>(p));
        if (mapValid_)
            cellMap_.add(particles_, particles_.size() - 1);
    }



    void map_particles() const
    {
        compact();
        cellMap_.add(particles_);
        mapValid_ = true;
    }
    void empty_map() { cellMap_.empty(); }


    NO_DISCARD std::size_t nbr_particles_in(box_t const& box) const
    {
        if (csrValid_)
        {
            std::size_t s = 0;
            forEachRow_(box, [&](int b, int e) { s += static_cast<std::size_t>(e - b); });
            for (auto i = csrPrefix_; i < particles_.size(); ++i)
                s += inBox_(particles_[i].iCell, box);
            return s;
        }
        ensureMap_();
        return cellMap_.size(box);
    }

    using cell_t = std::array<int, dim>;
    std::size_t nbr_particles_in(cell_t const& cell) const
    {
        if (csrValid_)
            return nbr_particles_in(box_t{cell, cell});
        ensureMap_();
        return cellMap_.size(cell);
    }

    /** nbr_particles_in(cell) for many cells of `box` (load balancer): sorted, the CSR row length
     * plus a histogram of the unsorted tail over `box` built once (O(cells + tail), not a tail
     * scan per cell); else the cell map. Same counts; valid while the array is not mutated.
     */
    class CellCounter
    {
    public:
        CellCounter(ParticleArray const& a, box_t const& box)
            : a_{a}
            , box_{box}
        {
            if (!a_.csrValid_)
            {
                a_.ensureMap_();
                return;
            }
            std::size_t n = 1;
            for (std::size_t d = 0; d < dim; ++d)
            {
                shape_[d] = box_.upper[d] - box_.lower[d] + 1;
                n *= static_cast<std::size_t>(shape_[d]);
            }
            tail_.assign(n, 0);
            for (auto i = a_.csrPrefix_; i < a_.particles_.size(); ++i)
                if (a_.inBox_(a_.particles_[i].iCell, box_))
                    ++tail_[key_(a_.particles_[i].iCell)];
        }

        std::size_t operator()(std::array<int, dim> const& cell) const
        {
            if (!a_.csrValid_)
                return a_.cellMap_.size(cell);
            std::size_t s = 0;
            if (a_.inBox_(cell, a_.box_))
            {
                auto const k = a_.key(cell);
                s            = static_cast<std::size_t>(a_.rowEnds()[k] - a_.csr_[k]);
            }
            if (a_.inBox_(cell, box_))
                return s + tail_[key_(cell)];
            for (auto i = a_.csrPrefix_; i < a_.particles_.size(); ++i)
                s += a_.particles_[i].iCell == cell;
            return s;
        }

    private:
        template<typename Cell>
        std::size_t key_(Cell const& cell) const
        {
            std::size_t k = 0;
            for (std::size_t d = 0; d < dim; ++d)
                k = k * static_cast<std::size_t>(shape_[d])
                    + static_cast<std::size_t>(cell[d] - box_.lower[d]);
            return k;
        }

        ParticleArray const& a_;
        box_t box_;
        std::array<int, dim> shape_{};
        std::vector<int> tail_;
    };

    NO_DISCARD CellCounter cellCounter(box_t const& box) const { return CellCounter{*this, box}; }

    void export_particles(box_t const& box, ParticleArray<dim>& dest) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles");
        if (csrValid_)
            return exportCSR_(box, [&](auto const& p) { dest.push_back(p); });
        ensureMap_();
        cellMap_.export_to(box, particles_, dest);
    }

    template<typename Fn>
    void export_particles(box_t const& box, ParticleArray<dim>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (Fn)");
        if (csrValid_)
            return exportCSR_(box, [&](auto const& p) { dest.push_back(fn(p)); });
        ensureMap_();
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (box, vector, Fn)");
        if (csrValid_)
            return exportCSR_(box, [&](auto const& p) { dest.push_back(fn(p)); });
        ensureMap_();
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Predicate>
    void export_particles(ParticleArray& dest, Predicate&& pred) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (Fn,vector)");
        if (csrValid_)
        {
            for_each_point(box_, [&](auto const& cell) {
                if (pred(cell))
                {
                    auto const k = key(cell);
                    for (int i = csr_[k]; i < rowEnds()[k]; ++i)
                        dest.push_back(particles_[i]);
                }
            });
            for (auto i = csrPrefix_; i < particles_.size(); ++i)
                if (inBox_(particles_[i].iCell, box_) and pred(particles_[i].iCell))
                    dest.push_back(particles_[i]);
            return;
        }
        ensureMap_();
        cellMap_.export_if(particles_.data(), dest, std::forward<Predicate>(pred));
    }


    template<typename Cell>
    void change_icell(Cell const& newCell, std::size_t particleIndex)
    {
        assert(!bagged_);
        csrValid_                       = false;
        auto oldCell                    = particles_[particleIndex].iCell;
        particles_[particleIndex].iCell = newCell;
        auto const box_is_valid         = box_.size() > 1;
        if (box_is_valid and mapValid_)
            cellMap_.update(particles_, particleIndex, oldCell);
    }

    // cell map side of change_icell for a particle whose iCell was already set: moves its index
    // from oldCell to newCell
    template<typename Cell>
    void map_icell_change(Cell const& oldCell, Cell const& newCell, std::size_t particleIndex)
    {
        csrValid_               = false;
        auto const box_is_valid = box_.size() > 1;
        if (box_is_valid and mapValid_)
            cellMap_.update(particleIndex, oldCell, newCell);
    }


    template<typename Predicate>
    auto partition(Predicate&& pred)
    {
        ensureMap_();
        csrValid_ = false;
        return cellMap_.partition(makeIndexRange(*this), std::forward<Predicate>(pred));
    }

    template<typename Range_t, typename Predicate>
    auto partition(Range_t&& range, Predicate&& pred)
    {
        ensureMap_();
        csrValid_      = false;
        auto const ret = cellMap_.partition(range, std::forward<Predicate>(pred));
        assert(ret.size() <= range.size());
        return ret;
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        ensureMap_();
        cellMap_.print(cell);
    }


    NO_DISCARD bool is_consistent() const
    {
        if (bagged_)
            return csrValid_
                   and checkBags([&](auto const& p) { return cellKey_(p.iCell); });
        if (!mapValid_)
            return csrValid_
                   and checkCSR([&](auto const& p) { return cellKey_(p.iCell); }, csrPrefix_);

        if (particles_.size() != cellMap_.size())
            return false;

        for (std::size_t pidx = 0; pidx < particles_.size(); ++pidx)
            if (!cellMap_(particles_[pidx].iCell).is_indexed(pidx))
                return false;

        return true;
    }

    void sortMapping() const
    {
        if (mapValid_)
            cellMap_.sort();
    }

    NO_DISCARD auto& vector() { return compact(), particles_; }
    NO_DISCARD auto& vector() const { return compact(), particles_; }

    auto& box() const { return box_; }


    // ---- sorted mode ----

    NO_DISCARD int nbrCells() const { return nCells_; }

    // the ion updater's kernel choice for this sorted array (IonUpdater::cellKernels_):
    // on = -1 not made yet, 0 the kernels walking the rows, 1 the cell kernels; correctors
    // counts the correctors since the array was created
    struct CellKernelChoice
    {
        std::int8_t on           = -1;
        std::uint32_t correctors = 0;
    };
    CellKernelChoice cellKernels;

    // row-major local flat index of a cell of box_ (last index fastest)
    template<typename Cell>
    NO_DISCARD int key(Cell const& cell) const
    {
        int k = 0;
        for (std::size_t d = 0; d < dim; ++d)
            k = k * shape_[d] + (cell[d] - box_.lower[d]);
        return k;
    }

    NO_DISCARD bool csrValid() const { return csrValid_; }
    // number of leading particles the CSR describes (0 if invalid)
    NO_DISCARD std::size_t csrPrefix() const { return csrValid_ ? csrPrefix_ : 0; }
    NO_DISCARD auto const& csr() const { return csr_; }
    // drops the CSR flag, keeping its content for sortBags(oldPrefix = csrPrefix())
    void invalidateCSR() { csrValid_ = false; }
    NO_DISCARD bool mapValid() const { return mapValid_; }
    // stop maintaining the cell map (rebuilt lazily by a query that needs it)
    void dropMap() { mapValid_ = false; }

    // ---- bags ----

    NO_DISCARD bool bagged() const { return bagged_; }
    // end of the live part of every row: row k = [csr()[k], rowEnds()[k]) (packed CSR: csr + 1)
    NO_DISCARD int const* rowEnds() const { return bagged_ ? ends_.data() : csr_.data() + 1; }

    NO_DISCARD std::size_t liveSize() const
    {
        if (!bagged_)
            return particles_.size();
        std::size_t s = particles_.size() - csrPrefix_;
        for (int k = 0; k < nCells_; ++k)
            s += static_cast<std::size_t>(ends_[k] - csr_[k]);
        return s;
    }

    /** remove the holes: rows moved forward in place, then the appends. A valid bag layout
     * becomes a valid packed CSR (rows keep their keys). No-op on a dense array.
     */
    void compact() const
    {
        if (bagged_)
            const_cast<This*>(this)->compact_();
    }

    /** bag sort. Requires the bag layout of before the in-place push (csr_, ends_ for
     * [0, oldPrefix)) and crossers = (index < oldPrefix, old key) of every particle whose cell
     * changed, ascending index; [oldPrefix, size) are unsorted appends. A crosser whose key
     * changed is copied out and its slot filled from the end of its row (holes stay at the row
     * tails); movers and appends go to a free slot of their row; a full row takes a slot from the
     * nearest row with one (cascade: one copy per non-empty row crossed). Key nCells() movers
     * (leavers) go to sc.leavers. Without a bag layout (first call, after a compaction) or when
     * no row within bagRadius has room: full layout with fresh slack. The array has no tail
     * afterwards (size() == csrPrefix()).
     */
    template<typename KeyFn>
    void sortBags(KeyFn&& keyOf, std::size_t const oldPrefix,
                  std::vector<std::pair<std::size_t, int>> const& crossers, SortScratch& sc)
    {
        sc.leavers.clear();
        if (!bagged_ or oldPrefix == 0)
        {
            layoutBags_(keyOf, sc);
            return;
        }
        auto const n = particles_.size();
        auto const K = nCells_;
        assert(static_cast<std::size_t>(csr_[K]) == oldPrefix);
        auto& buf  = sc.buf;
        auto& bufK = sc.bufKeys;
        auto& vac  = sc.vacated;
        buf.clear();
        bufK.clear();
        vac.clear();

        // 1. copy the movers out (ascending index = ascending row)
        constexpr std::size_t pf = 16;
        auto const nCross        = crossers.size();
        for (std::size_t c = 0; c < nCross; ++c)
        {
            if (c + pf < nCross)
                __builtin_prefetch(&particles_[crossers[c + pf].first]);
            auto const& [idx, ok] = crossers[c];
            int const k           = keyOf(particles_[idx]);
            if (k == ok)
                continue;
            buf.push_back(particles_[idx]);
            bufK.push_back(k);
            vac.push_back(static_cast<int>(idx));
        }
        // 2. close the vacated slots of each row with the row's last live particles
        int row            = 0;
        for (std::size_t a = 0; a < vac.size();)
        {
            while (csr_[row + 1] <= vac[a])
                ++row;
            std::size_t b = a;
            while (b < vac.size() and vac[b] < csr_[row + 1])
                ++b;
            // vac[a, b) in row, ascending
            int const e  = ends_[row];
            int const ne = e - static_cast<int>(b - a);
            std::size_t hi = b; // vacated slots >= ne: vac[firstHi, b)
            while (hi > a and vac[hi - 1] >= ne)
                --hi;
            std::size_t h = a, sk = hi;
            for (int d = ne; d < e and h < hi; ++d)
            {
                if (sk < b and vac[sk] == d)
                {
                    ++sk;
                    continue;
                }
                particles_[vac[h++]] = particles_[d];
            }
            ends_[row] = ne;
            a          = b;
        }
        // 3. the appends
        for (auto i = oldPrefix; i < n; ++i)
        {
            buf.push_back(particles_[i]);
            bufK.push_back(keyOf(particles_[i]));
        }
        // 4. place
        auto const m = buf.size();
        for (std::size_t j = 0; j < m; ++j)
        {
            if (j + pf < m and bufK[j + pf] < K)
                __builtin_prefetch(&particles_[std::min<std::size_t>(
                                       static_cast<std::size_t>(ends_[bufK[j + pf]]), n - 1)],
                                   1);
            int const k = bufK[j];
            if (k == K)
            {
                sc.leavers.push_back(buf[j]);
                continue;
            }
            if (ends_[k] == csr_[k + 1] and !makeRoom_(k))
            {
                // no free slot near row k: the rest goes back as appends, full layout
                particles_.resize(static_cast<std::size_t>(csr_[K]));
                csr_[K + 1] = csr_[K];
                csrPrefix_  = static_cast<std::size_t>(csr_[K]);
                for (auto jj = j; jj < m; ++jj)
                    if (bufK[jj] == K)
                        sc.leavers.push_back(buf[jj]);
                    else
                        particles_.push_back(buf[jj]);
                layoutBags_(keyOf, sc, /*keepLeavers=*/true);
                return;
            }
            particles_[static_cast<std::size_t>(ends_[k]++)] = buf[j];
        }
        particles_.resize(static_cast<std::size_t>(csr_[K]));
        finishBags_();
    }

    // every row [csr_[k], ends_[k]) holds key k only, capacities nest, no tail
    template<typename KeyFn>
    NO_DISCARD bool checkBags(KeyFn&& keyOf) const
    {
        auto const K = nCells_;
        if (!bagged_ or csr_.size() != static_cast<std::size_t>(K + 2) or ends_.size() != static_cast<std::size_t>(K)
            or csr_[0] != 0 or csr_[K] != csr_[K + 1]
            or static_cast<std::size_t>(csr_[K]) != csrPrefix_ or csrPrefix_ > particles_.size())
            return false;
        for (int k = 0; k < K; ++k)
        {
            if (csr_[k] > ends_[k] or ends_[k] > csr_[k + 1])
                return false;
            for (int i = csr_[k]; i < ends_[k]; ++i)
                if (keyOf(particles_[i]) != k)
                    return false;
        }
        return true;
    }

    // every particle of [0, prefix) has keyOf == its row, and csr_ ends at prefix
    template<typename KeyFn>
    NO_DISCARD bool checkCSR(KeyFn&& keyOf, std::size_t const prefix) const
    {
        if (csr_.size() != static_cast<std::size_t>(nCells_ + 2)
            or static_cast<std::size_t>(csr_[nCells_ + 1]) != prefix or csr_[0] != 0)
            return false;
        for (int k = 0; k <= nCells_; ++k)
            for (int i = csr_[k]; i < csr_[k + 1]; ++i)
                if (keyOf(particles_[i]) != k)
                    return false;
        return true;
    }


    /** fn(first, last) over the dense runs of live particles, in storage order: the rows
     * (consecutive rows without a hole between them make one run), then the appends. The order
     * compact() would leave, without moving anything.
     */
    template<typename Fn>
    void forEachLiveRun(Fn&& fn) const
    {
        std::size_t s = 0, e = 0; // current run [s, e)
        auto const add = [&](std::size_t const b, std::size_t const en) {
            if (b == en)
                return;
            if (b != e)
            {
                if (e > s)
                    fn(s, e);
                s = b;
            }
            e = en;
        };
        if (bagged_)
            for (int k = 0; k < nCells_; ++k)
                add(static_cast<std::size_t>(csr_[k]), static_cast<std::size_t>(ends_[k]));
        add(bagged_ ? csrPrefix_ : 0, particles_.size());
        if (e > s)
            fn(s, e);
    }

    /** packed CSR layout (rows without slack) of every live particle, in their current order
     * within a row; key nCells() particles stay after the rows as unsorted appends. For an
     * array that reaches the push unsorted (new, restarted, regridded): a packed CSR is also
     * what a compaction leaves, so a restarted run lays its arrays out exactly as the
     * uninterrupted run whose restart dump compacted them.
     */
    template<typename KeyFn>
    void layoutRows(KeyFn&& keyOf, SortScratch& sc)
    {
        layoutBags_(keyOf, sc, /*keepLeavers=*/false, /*slots=*/0);
        for (auto const& p : sc.leavers)
            particles_.push_back(p);
        sc.leavers.clear();
    }


private:
    // full layout of every live particle (rows, then the appends): out-of-place counting
    // sort into rows of capacity count + slots (bags, or a packed CSR with slots 0), leavers
    // (key nCells_) to sc.leavers
    template<typename KeyFn>
    void layoutBags_(KeyFn&& keyOf, SortScratch& sc, bool const keepLeavers = false,
                     int const slots = bagSlots)
    {
        if (!keepLeavers)
            sc.leavers.clear();
        auto const K = nCells_;
        auto& cnt    = sc.counts;
        cnt.assign(K + 1, 0);
        auto& keys = sc.keys;
        keys.clear();
        forEachLive_([&](std::size_t i) {
            int const k = keyOf(particles_[i]);
            keys.push_back(k);
            ++cnt[k];
        });
        // new layout in scratch: the second pass walks the old one
        auto& start = sc.csr;
        auto& end   = sc.ends;
        start.resize(K + 2);
        end.resize(K);
        start[0] = 0;
        for (int k = 0; k < K; ++k)
        {
            int const c  = cnt[k];
            start[k + 1] = start[k] + c + (c > 0 ? slots : 0);
            end[k]       = start[k];
        }
        start[K + 1] = start[K];
        auto& dst    = sc.particles;
        dst.resize(static_cast<std::size_t>(start[K]));
        std::size_t j = 0;
        forEachLive_([&](std::size_t i) {
            int const k = keys[j++];
            if (k == K)
                sc.leavers.push_back(particles_[i]);
            else
                dst[static_cast<std::size_t>(end[k]++)] = particles_[i];
        });
        particles_.swap(dst);
        csr_.swap(start);
        ends_.swap(end);
        bagged_ = slots > 0;
        finishBags_();
    }

    void finishBags_()
    {
        auto const K = nCells_;
        csr_[K + 1]  = csr_[K];
        csrPrefix_   = static_cast<std::size_t>(csr_[K]);
        assert(particles_.size() == csrPrefix_);
        csrValid_ = true;
        mapValid_ = false;
    }

    // fn(i) for every live slot: rows (bagged) or [0, csrPrefix_), then [csrPrefix_, size)
    template<typename Fn>
    void forEachLive_(Fn&& fn) const
    {
        std::size_t from = 0;
        if (bagged_)
        {
            for (int k = 0; k < nCells_; ++k)
                for (int i = csr_[k]; i < ends_[k]; ++i)
                    fn(static_cast<std::size_t>(i));
            from = csrPrefix_;
        }
        for (auto i = from; i < particles_.size(); ++i)
            fn(i);
    }

    /** one free slot at the end of the full row k, taken from the nearest row j with a free
     * slot (|j - k| <= bagRadius); the rows between shift by one slot toward j (one copy per
     * non-empty row). Past the last row the capacity end grows (storage resized). false: none.
     */
    bool makeRoom_(int const k)
    {
        auto const K = nCells_;
        int const R  = bagRadius;
        int right = -1, left = -1;
        for (int j = k + 1; j <= std::min(K, k + R); ++j)
            if (j == K or ends_[j] < csr_[j + 1])
            {
                right = j;
                break;
            }
        for (int j = k - 1; j >= std::max(0, k - R); --j)
            if (ends_[j] < csr_[j + 1])
            {
                left = j;
                break;
            }
        if (right < 0 and left < 0)
            return false;
        if (right == K and left >= 0)
            right = -1; // reuse a hole rather than grow the storage
        if (right >= 0 and (left < 0 or right - k <= k - left))
        {
            if (right == K)
            {
                // grow the capacity end by one slot
                ++csr_[K];
                ++csr_[K + 1];
                if (particles_.size() < static_cast<std::size_t>(csr_[K]))
                    particles_.resize(static_cast<std::size_t>(csr_[K]));
            }
            else
            {
                // row right: first live particle to its free end slot
                if (ends_[right] > csr_[right])
                {
                    particles_[ends_[right]] = particles_[csr_[right]];
                }
                ++ends_[right];
                ++csr_[right];
            }
            for (int m = (right == K ? K : right) - 1; m > k; --m)
            {
                // full row m: first to the slot freed after its end
                if (ends_[m] > csr_[m])
                {
                    particles_[ends_[m]] = particles_[csr_[m]];
                }
                ++ends_[m];
                ++csr_[m];
            }
            return true;
        }
        // left: rows left+1 .. k shift down one slot (last live particle to the slot before)
        for (int m = left + 1; m <= k; ++m)
        {
            if (ends_[m] > csr_[m])
            {
                particles_[csr_[m] - 1] = particles_[ends_[m] - 1];
            }
            --ends_[m];
            --csr_[m];
        }
        return true;
    }

    void compact_()
    {
        auto const K   = nCells_;
        std::size_t w  = 0;
        auto const old = particles_.size();
        for (int k = 0; k < K; ++k)
        {
            int const s = csr_[k], e = ends_[k];
            csr_[k]     = static_cast<int>(w);
            for (int i = s; i < e; ++i)
                particles_[w++] = particles_[static_cast<std::size_t>(i)];
        }
        auto const rowsEnd = w;
        for (auto i = csrPrefix_; i < old; ++i)
            particles_[w++] = particles_[i];
        particles_.resize(w);
        csr_[K]     = static_cast<int>(rowsEnd);
        csr_[K + 1] = static_cast<int>(rowsEnd);
        csrPrefix_  = rowsEnd;
        bagged_     = false;
        // csrValid_ unchanged: a valid bag layout is now a valid packed CSR
    }

    void initShape_()
    {
        nCells_ = 1;
        for (std::size_t d = 0; d < dim; ++d)
        {
            shape_[d] = box_.upper[d] - box_.lower[d] + 1;
            nCells_ *= shape_[d];
        }
    }

    template<typename Cell>
    bool inBox_(Cell const& cell, box_t const& box) const
    {
        for (std::size_t d = 0; d < dim; ++d)
            if (cell[d] < box.lower[d] or cell[d] > box.upper[d])
                return false;
        return true;
    }

    template<typename Cell>
    int cellKey_(Cell const& cell) const
    {
        return inBox_(cell, box_) ? key(cell) : nCells_;
    }

    void ensureMap_() const
    {
        if (mapValid_)
            return;
        compact();
        cellMap_.empty();
        cellMap_.add(particles_);
        mapValid_ = true;
    }

    // fn(rowBegin, rowEnd) for each contiguous key run of box clipped to box_
    template<typename Fn>
    void forEachRow_(box_t const& query, Fn&& fn) const
    {
        std::array<int, dim> lo, hi;
        for (std::size_t d = 0; d < dim; ++d)
        {
            lo[d] = std::max(query.lower[d], box_.lower[d]);
            hi[d] = std::min(query.upper[d], box_.upper[d]);
            if (lo[d] > hi[d])
                return;
        }
        constexpr std::size_t last = dim - 1;
        auto row                   = [&](std::array<int, dim> c) {
            c[last]      = lo[last];
            int const k0 = key(c);
            int const k1 = k0 + (hi[last] - lo[last]);
            if (bagged_)
                for (int k = k0; k <= k1; ++k)
                    fn(csr_[k], ends_[k]);
            else
                fn(csr_[k0], csr_[k1 + 1]);
        };
        if constexpr (dim == 1)
            row(lo);
        else if constexpr (dim == 2)
        {
            for (int i = lo[0]; i <= hi[0]; ++i)
                row({i, lo[1]});
        }
        else
        {
            for (int i = lo[0]; i <= hi[0]; ++i)
                for (int j = lo[1]; j <= hi[1]; ++j)
                    row({i, j, lo[2]});
        }
    }

    template<typename Push>
    void exportCSR_(box_t const& box, Push&& push) const
    {
        forEachRow_(box, [&](int b, int e) {
            for (int i = b; i < e; ++i)
                push(particles_[i]);
        });
        for (auto i = csrPrefix_; i < particles_.size(); ++i)
            if (inBox_(particles_[i].iCell, box))
                push(particles_[i]);
    }

    Vector particles_;
    box_t box_;
    mutable CellMap_t cellMap_;

    std::array<int, dim> shape_{};
    int nCells_ = 0;
    std::vector<int> csr_;
    std::size_t csrPrefix_ = 0;
    bool csrValid_         = false;
    mutable bool mapValid_ = true;
    // bags: row k = live [csr_[k], ends_[k]), holes [ends_[k], csr_[k + 1]); csr_[nCells_] =
    // csr_[nCells_ + 1] = csrPrefix_ = end of the rows' capacity; [csrPrefix_, size) = appends
    bool bagged_ = false;
    std::vector<int> ends_;
};

} // namespace PHARE::core


namespace PHARE
{
namespace core
{


    template<std::size_t dim, bool OwnedState = true>
    struct ContiguousParticles
    {
        static constexpr bool is_contiguous    = true;
        static constexpr std::size_t dimension = dim;
        using ContiguousParticles_             = ContiguousParticles<dim, OwnedState>;

        template<typename T>
        using container_t = std::conditional_t<OwnedState, std::vector<T>, Span<T>>;

        template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
        ContiguousParticles(std::size_t s)
            : iCell(s * dim)
            , delta(s * dim)
            , weight(s)
            , v(s * 3)
        {
        }

        template<typename Container_int, typename Container_double>
        ContiguousParticles(Container_int&& _iCell, Container_double&& _delta,
                            Container_double&& _weight, Container_double&& _v)
            : iCell{_iCell}
            , delta{_delta}
            , weight{_weight}
            , v{_v}
        {
        }

        NO_DISCARD std::size_t size() const { return weight.size(); }

        template<std::size_t S, typename T>
        NO_DISCARD static std::array<T, S>* _array_cast(T const* array)
        {
            return reinterpret_cast<std::array<T, S>*>(const_cast<T*>(array));
        }

        template<typename Return>
        NO_DISCARD Return _to(std::size_t i)
        {
            return {
                *const_cast<double*>(weight.data() + i),     //
                *_array_cast<dim>(iCell.data() + (dim * i)), //
                *_array_cast<dim>(delta.data() + (dim * i)), //
                *_array_cast<3>(v.data() + (3 * i)),
            };
        }

        NO_DISCARD auto copy(std::size_t i) { return _to<Particle<dim>>(i); }
        NO_DISCARD auto view(std::size_t i) { return _to<ParticleView<dim>>(i); }

        NO_DISCARD auto operator[](std::size_t i) const { return view(i); }
        NO_DISCARD auto operator[](std::size_t i) { return view(i); }

        struct iterator
        {
            iterator(ContiguousParticles_* particles)
            {
                for (std::size_t i = 0; i < particles->size(); i++)
                    views.emplace_back((*particles)[i]);
            }

            iterator& operator++()
            {
                ++curr_pos;
                return *this;
            }

            NO_DISCARD bool operator!=(iterator const& other) const
            {
                return curr_pos != views.size();
            }
            NO_DISCARD auto& operator*() { return views[curr_pos]; }
            NO_DISCARD auto& operator*() const { return views[curr_pos]; }

            std::size_t curr_pos = 0;
            std::vector<ParticleView<dim>> views;
        };

        NO_DISCARD auto as_tuple() { return std::forward_as_tuple(weight, iCell, delta, v); }
        NO_DISCARD auto as_tuple() const { return std::forward_as_tuple(weight, iCell, delta, v); }

        NO_DISCARD auto begin() { return iterator(this); }
        NO_DISCARD auto cbegin() const { return iterator(this); }

        NO_DISCARD auto end() { return iterator(this); }
        NO_DISCARD auto cend() const { return iterator(this); }

        container_t<int> iCell;
        container_t<double> delta;
        container_t<double> weight, v;
    };


    template<std::size_t dim>
    using ContiguousParticlesView = ContiguousParticles<dim, /*OwnedState=*/false>;

} // namespace core
} // namespace PHARE


#endif
