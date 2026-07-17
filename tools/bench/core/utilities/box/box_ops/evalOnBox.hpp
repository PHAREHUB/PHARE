#define OMPI_SKIP_MPICXX
#include <cassert>
#include <tuple>
#include <vector>
#include <array>
#include <chrono>
// #include <span>

#include "mkn/kul/log.hpp"
#include "mkn/kul/time.hpp"
#include "mkn/avx/vector.hpp"

#include "phare_core.hpp"
#include "core/utilities/span.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"


constexpr static std::size_t dim    = 3;
constexpr static std::size_t interp = 1;

constexpr static std::size_t C = 16;
constexpr static std::size_t N = 5000;

using T            = double;
using PHARE_Types  = PHARE::core::PHARE_Types<dim, 1>;
using GridLayout_t = TestGridLayout<typename PHARE_Types::GridLayout_t>;
using Array_t      = PHARE::core::NdArrayVector<dim, T>;
using Grid_t       = PHARE::core::Grid<Array_t, PHARE::core::HybridQuantity::Scalar>;


namespace PHARE::core
{

template<typename Array_t>
struct BoxRow
{
    Array_t& arr;
    Box<std::uint32_t, dim> box;
    std::uint32_t k;
    std::uint32_t j = box.lower[1];
    std::uint32_t s = box.upper[2] - box.lower[2] + 1;
    Span<double> row{&arr(k, j, box.lower[2]), s};

    BoxRow& operator++()
    {
        row = Span<double>{&arr(k, ++j, box.lower[2]), s};
        return *this;
    }
    auto& operator*() { return row; }

    bool operator==(BoxRow const& that) const { return j == that.j; }
    bool operator!=(BoxRow const& that) const { return j != that.j; }

    double const* const last_domain_p_1 = &arr(box.upper);

    void next()
    {
        j   = box.lower[1];
        row = Span<double>{&arr(++k, j, box.lower[2]), s};
    }
};

template<typename Array_t>
struct BoxSlab
{
    using BoxRow_t = BoxRow<Array_t>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    bool _end = false;
    BoxRow<Array_t> br{arr, box, _end ? box.upper[0] : box.lower[0]};

    BoxRow_t begin() const
    {
        assert(!_end);
        return br;
    }
    BoxRow_t end() const { return {arr, box, box.upper[0], box.upper[1] + 1, 0, {0, 0}}; }

    void next() { br.next(); }
};

template<typename Array_t>
struct BoxSlabber
{
    using BoxSlab_t = BoxSlab<Array_t>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    bool _end       = false;
    std::uint32_t k = _end ? box.upper[0] + 1 : box.lower[0];

    BoxSlab_t slab{arr, box, _end};

    bool operator==(BoxSlabber const& that) const { return k == that.k; }
    bool operator!=(BoxSlabber const& that) const { return k != that.k; }

    BoxSlabber& operator++()
    {
        slab.next();
        ++k;
        return *this;
    }
    auto& operator*() { return slab; }
};

template<typename Array_t>
struct BoxSpan
{
    using BoxSlabber_t = BoxSlabber<Array_t>;

    Box<std::uint32_t, dim> box;
    Array_t& arr;

    BoxSlabber_t b{arr, box};
    BoxSlabber_t e{arr, box, true};

    auto begin() { return b; }
    auto end() { return e; }
    auto begin() const { return b; }
    auto end() const { return e; }
};


template<typename Array_t>
auto make_box_span(Box<std::uint32_t, dim> const box, Array_t& arr)
{
    return BoxSpan{box, arr};
}

template<typename Array_t>
auto make_box_span(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    return BoxSpan{box, arr};
}


std::uint64_t static now()
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

struct Timer
{
    // Timer() {}
    ~Timer()
    {
        auto const total = now() - s;
        PHARE_LOG_LINE_SS("RUN: " << total << " ns");
        PHARE_LOG_LINE_SS("AVG: " << (total / div) << " ns");
    }

    std::size_t div = 1;
    std::size_t const s{now()};
};


} // namespace PHARE::core
