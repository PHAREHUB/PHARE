#ifndef PHARE_CORE_DEBUG_HPP
#define PHARE_CORE_DEBUG_HPP

/*
    For debugging field values and overlaps during execution
      If building in Release mode, you should  use the compile flag
        -DPHARE_FORCE_LOG_LINE=1
      To enable logging, while keeping asserts disabled

    This is only expected to work for serial operations!
*/


#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"

#include <cassert>
#include <sstream>
#include <unordered_map>

#if !defined(PHARE_DEBUGGERINO)
#define PHARE_DEBUGGERINO 0
#endif


#if PHARE_DEBUGGERINO
// AVAILABLE FUNCTIONS

// clang-format off

// Tag scope to make function path like /main/func1/func1/func3/etc
#define PHARE_DEBUG_SCOPE(key)  \
    PHARE::core::debug_scope _debug_scope_{key}

// Set current level data
#define PHARE_DEBUG_SET(time, lvlnum, hierarchy) \
    PHARE::core::Debuggerino::INSTANCE().set(time, lvlnum, hierarchy)

// Monitor FieldBox operations for particular fields and particular amr boxes
//  uses the debugging map up the top of this file
#define PHARE_DEBUG_FIELDS(dst, src)  \
    PHARE::core::debug_if_active_fields(dst, src)

// Monitor Field values for particular fields and particular amr boxes
//  uses the debugging map up the top of this file
#define PHARE_DEBUG_FIELD(field, layout)  \
    PHARE::core::debug_if_active_field(field, layout)

// Monitor Field values for all fields and particular amr boxes
//  uses model views
#define PHARE_DEBUG_ALL_FIELDS(views)  \
    PHARE::core::debug_all_fields(views)

// Check overlaps from model views
#define PHARE_DEBUG_CHECK_ALL_OVERLAPS(views)  \
    PHARE::core::check_all_fields(views)

// Check overlaps from rm and SAMRAI::PatchLevel
#define PHARE_DEBUG_CHECK_LEVEL(GridLayoutT, rm, level) \
    PHARE::core::check_level<GridLayoutT>(rm, level)

#endif
// clang-format on

namespace PHARE::core
{

std::size_t constexpr inline static DEBUG_DIM = 2;

struct ActiveFieldInfo
{
    double time;
    Box<int, DEBUG_DIM> box;
    std::string dst_name, src_name = dst_name;
};

// some examples provided
//  To see all available paths. export env var `PHARE_DEBUG_PRINT_ALL_FIELDS=1`
std::unordered_map<std::string, std::vector<ActiveFieldInfo>> static inline const debugging{
    {
        "Simulator/advance/level/1/SolverPPC/predictor1/",
        {
            // {.225, Box<int, DEBUG_DIM>{{62, 36}, {64, 42}}, "rho"},
            // {.225, Box<int, DEBUG_DIM>{{62, 36}, {64, 42}}, "EM_B_x"},
        } //
    },
    {
        "Simulator/advance/HyHyMessStrat/after_sync/",
        {
            // {.225, Box<int, DEBUG_DIM>{{62, 36}, {64, 42}}, "EM_B_x"},
            // {.225, Box<int, DEBUG_DIM>{{62, 36}, {64, 42}}, "EM_B_y"},
            // {.225, Box<int, DEBUG_DIM>{{62, 36}, {64, 42}}, "EM_B_z"},
        } //
    } //
};

// If empty checks ALL FIELDS!!
std::vector<std::string> static inline const active_overlap_fields{{"EM_B_x"}};


struct debug_scope;

struct Debuggerino
{
    static Debuggerino& INSTANCE();

    Debuggerino() = default;

    void set(double const t, int l, auto h)
    {
        time              = t;
        level             = l;
        auto domainBoxVec = h->domainBox();
        domainBox         = {Point{ConstArray<int, DEBUG_DIM>(0)},
                             Point{for_N<DEBUG_DIM, for_N_R_mode::make_array>(
                         [&](auto i) { return domainBoxVec[i]; })}};
        if (l)
            domainBox.upper = (domainBox.upper + 1) * (l + 1) - 1;
    }

    double time = 0;
    int level   = 0;
    Box<int, DEBUG_DIM> domainBox;
    debug_scope* stack_ptr = nullptr;
};

struct Shifter
{
    template<typename T, std::size_t D>
    auto make_shift_for(PHARE::core::Box<T, D> const& box) const
    {
        auto const span = domainBox.upper + 1;
        auto const mid  = for_N<D, for_N_R_mode::make_array>([&](auto i) { return span[i] / 2; });

        if constexpr (D == 1)
        {
            auto shifts = PHARE::core::for_N<1, PHARE::core::for_N_R_mode::make_array>(
                [&](auto i) { return PHARE::core::Point<int, D>{0}; });
            PHARE::core::for_N<D>([&](auto i) {
                int const shift = box.upper[i] < mid[i] ? 1 : -1;
                shifts[i][i]    = span[i] * shift;
            });

            return shifts;
        }
        if constexpr (D == 2)
        {
            auto shifts = PHARE::core::for_N<3, PHARE::core::for_N_R_mode::make_array>(
                [&](auto i) { return PHARE::core::Point<int, D>{0, 0}; });
            PHARE::core::for_N<D>([&](auto i) {
                int const shift = box.upper[i] < mid[i] ? 1 : -1;
                shifts[i][i]    = span[i] * shift;
            });

            shifts[2] = {shifts[0][0], shifts[1][1]};

            return shifts;
        }
        if constexpr (D == 3)
        {
            auto shifts = PHARE::core::for_N<7, PHARE::core::for_N_R_mode::make_array>(
                [&](auto i) { return PHARE::core::Point<int, D>{0, 0, 0}; });
            PHARE::core::for_N<D>([&](auto i) {
                int const shift = box.upper[i] < mid[i] ? 1 : -1;
                shifts[i][i]    = span[i] * shift;
            });

            shifts[3] = {shifts[0][0], shifts[1][1], 0};
            shifts[4] = {0, shifts[1][1], shifts[2][2]};
            shifts[5] = {shifts[0][0], 0, shifts[2][2]};
            shifts[6] = {shifts[0][0], shifts[1][1], shifts[2][2]};

            return shifts;
        }
    }


    Box<int, DEBUG_DIM> domainBox;
};


struct debug_scope
{
    debug_scope(std::string const& key);

    ~debug_scope();

    static debug_scope const* root_parent_from(auto const& self)
    {
        if (self.parent)
            return root_parent_from(*self.parent);
        return &self;
    }


    auto full_path() const
    {
        auto const self = this;
        auto iter       = root_parent_from(*this);
        std::stringstream ss;

        ss << iter->key;

        while (iter and iter != this)
            if ((iter = iter->child ? iter->child : nullptr))
                ss << iter->key;

        return ss.str();
    }

    std::string key;
    debug_scope* parent = nullptr;
    debug_scope* child  = nullptr;
};


void debug_fields(auto& dst, auto& src, auto& k, auto& v)
{
    auto& debugger           = Debuggerino::INSTANCE();
    auto& scope              = *debugger.stack_ptr;
    auto const path          = scope.full_path();
    auto const dst_ghost_box = shift(dst.amr_ghost_box, src.offset_ * -1);
    auto const src_ghost_box = src.amr_ghost_box;

    PHARE_LOG_LINE_SS("debugging: " << path << " / time: " << v.time
                                    << " / field: " << dst.field.name());

    auto src_it = src.lcl_box.begin();
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++src_it, ++dst_it)
    {
        auto const amr_idx = (*dst_it + dst_ghost_box.lower).as_signed();
        if (isIn(amr_idx, v.box))
        {
            auto const dv = to_string_with_precision(dst.field(*dst_it), 17);
            auto const sv = to_string_with_precision(src.field(*src_it), 17);

            std::cout << "amr_idx=" << amr_idx << " dst(" << dv << ") = src(" << sv << ")"
                      << std::endl;

            if (auto src_amr_idx = (*src_it + src_ghost_box.lower).as_signed();
                amr_idx != src_amr_idx)
            {
                PHARE_LOG_LINE_SS("amr_idx MISMATCH dst=" << amr_idx << " src=" << src_amr_idx);
            }
        }
    }
}

void debug_if_active_fields(auto const& dst, auto const& src)
{
    if constexpr (DEBUG_DIM == std::decay_t<decltype(dst)>::dimension)
    {
        assert(Debuggerino::INSTANCE().stack_ptr);
        auto& debugger           = Debuggerino::INSTANCE();
        debug_scope& scope       = *debugger.stack_ptr;
        auto const path          = scope.full_path();
        auto const dst_ghost_box = shift(dst.amr_ghost_box, src.offset_ * -1);
        auto const dst_amr_box   = shift(
            Box<int, DEBUG_DIM>{dst.lcl_box.lower.as_signed(), dst.lcl_box.upper.as_signed()},
            dst_ghost_box.lower);

        if (get_env("PHARE_DEBUG_PRINT_ALL_FIELDS"))
        {
            PHARE_LOG_LINE_SS(path << " / time: " << debugger.time
                                   << " / field: " << dst.field.name() << " " << dst_amr_box);
        }

        for (auto [k, vs] : debugging)
            for (auto const& v : vs)
                if (float_equals(debugger.time, v.time) and path == k
                    and v.dst_name == dst.field.name() and v.box * dst_amr_box)
                    debug_fields(dst, src, k, v);
    }
}



void debug_field(auto& field, auto& layout, auto& field_info)
{
    auto& debugger           = Debuggerino::INSTANCE();
    auto& scope              = *debugger.stack_ptr;
    auto const path          = scope.full_path();
    auto const dst_ghost_box = layout.AMRGhostBoxFor(field.physicalQuantity());

    PHARE_LOG_LINE_SS("debugging: " << path << " / time: " << field_info.time
                                    << " / field: " << field.name() << " patch" << layout.AMRBox());

    if (auto const overlap = field_info.box * dst_ghost_box)
        for (auto const& bix : layout.AMRToLocal(*overlap))
        {
            auto const amr_idx = (bix + dst_ghost_box.lower).as_signed();
            if (isIn(amr_idx, field_info.box))
            {
                auto const dv = to_string_with_precision(field(*bix), 17);
                std::cout << "amr_idx=" << amr_idx << " dst(" << dv << ") " << std::endl;
            }
        }
}

void debug_if_active_field(auto const& field, auto const& layout)
{
    if constexpr (DEBUG_DIM == std::decay_t<decltype(field)>::dimension)
    {
        assert(Debuggerino::INSTANCE().stack_ptr);
        auto& debugger           = Debuggerino::INSTANCE();
        debug_scope& scope       = *debugger.stack_ptr;
        auto const path          = scope.full_path();
        auto const dst_ghost_box = layout.AMRGhostBoxFor(field.physicalQuantity());

        if (get_env("PHARE_DEBUG_PRINT_ALL_FIELDS"))
        {
            PHARE_LOG_LINE_SS(path << " / time: " << debugger.time << " / field: " << field.name()
                                   << " " << dst_ghost_box);
        }

        for (auto [stack_path, field_infos] : debugging)
            for (auto const& field_info : field_infos)
                if (float_equals(debugger.time, field_info.time) and path == stack_path
                    and field_info.dst_name == field.name() and field_info.box * dst_ghost_box)
                    debug_field(field, layout, field_info);
    }
}


void debug_all_fields(auto& views)
{
    using ResourcesManager_t = std::decay_t<decltype(*views.model().resourcesManager)>;
    using FieldData_t        = typename ResourcesManager_t::UserField_t::patch_data_type;
    auto const& rm           = *views.model().resourcesManager;

    for (auto [name, info] : rm.all_resources())
    {
        for (auto& state : views)
        {
            auto data = state.patch->getPatchData(info.id);
            if (auto field_ptr = dynamic_cast<FieldData_t*>(data.get()))
                debug_if_active_field(*field_ptr->getPointer(), state.layout);
        }
    }
}

void check_fields(auto const& f0, auto const& f1)
{
    auto const name = f0->field.name();

    bool const active
        = active_overlap_fields.size() == 0
          || any(active_overlap_fields, [&](auto const& field_name) { return name == field_name; });

    if (not active)
        return;

    auto& debugger     = Debuggerino::INSTANCE();
    auto& scope        = *debugger.stack_ptr;
    auto const path    = scope.full_path();
    auto const layout0 = f0->gridLayout;
    auto const layout1 = f1->gridLayout;
    auto const pq      = f0->field.physicalQuantity();
    auto const gb0     = layout0.AMRGhostBoxFor(pq);
    auto const gb1     = layout1.AMRGhostBoxFor(pq);
    auto const shifter = Shifter{debugger.domainBox};


    // periodic checks
    for (auto const& offset : shifter.make_shift_for(layout0.AMRBox()))
    {
        auto const ghostbox = shift(gb0, offset);

        if (auto const overlap_ = ghostbox * gb1)
        {
            auto const overlap = *overlap_;
            auto const lb0     = layout0.AMRToLocal(shift(overlap, offset * -1));
            auto const lb1     = layout1.AMRToLocal(overlap);
            auto itamr         = overlap.begin();
            auto it0           = lb0.begin();
            auto it1           = lb1.begin();
            for (; it0 != lb0.end(); ++it0, ++it1, ++itamr)
            {
                auto const v0 = f0->field(*it0);
                auto const v1 = f1->field(*it1);
                if (!float_equals(v0, v1))
                {
                    PHARE_LOG_LINE_SS("periodic fail " << path << " " << f0->field.name()
                                                       << " time " << debugger.time);
                    PHARE_LOG_LINE_SS(layout0.AMRBox() << " " << layout1.AMRBox());
                    PHARE_LOG_LINE_SS((*itamr) << " " << to_string_with_precision(v0, 17) << " "
                                               << to_string_with_precision(v1, 17));
                }
            }
        }
    }

    // non-periodic check
    if (auto overlap_ = gb0 * gb1)
    {
        auto const overlap = *overlap_;
        auto const lb0     = layout0.AMRToLocal(overlap);
        auto const lb1     = layout1.AMRToLocal(overlap);

        auto itamr = overlap.begin();
        auto it0   = lb0.begin();
        auto it1   = lb1.begin();
        for (; it0 != lb0.end(); ++it0, ++it1, ++itamr)
        {
            auto const v0 = f0->field(*it0);
            auto const v1 = f1->field(*it1);
            if (!float_equals(v0, v1))
            {
                PHARE_LOG_LINE_SS("fail " << path << " " << f0->field.name() << " time "
                                          << debugger.time);
                PHARE_LOG_LINE_SS(layout0.AMRBox() << " " << layout1.AMRBox());
                PHARE_LOG_LINE_SS((*itamr) << " " << to_string_with_precision(v0, 17) << " "
                                           << to_string_with_precision(v1, 17));
            }
        }
    }
}

void check_all_fields(auto& views)
{
    using ResourcesManager_t = std::decay_t<decltype(*views.model().resourcesManager)>;
    using FieldData_t        = typename ResourcesManager_t::UserField_t::patch_data_type;

    if constexpr (DEBUG_DIM == FieldData_t::dimension)
    {
        auto const& rm     = *views.model().resourcesManager;
        auto& debugger     = Debuggerino::INSTANCE();
        auto const& states = views.states;

        for (auto [name, info] : rm.all_resources())
            for (std::size_t i = 0; i < states.size() - 1; ++i)
                for (std::size_t j = i + 1; j < states.size(); ++j)
                {
                    auto const& p0 = *states[i].patch;
                    auto const& p1 = *states[j].patch;
                    auto const& d0 = p0.getPatchData(info.id);
                    auto const& d1 = p1.getPatchData(info.id);
                    auto f0        = dynamic_cast<FieldData_t*>(d0.get());
                    if (f0)
                        check_fields(f0, dynamic_cast<FieldData_t*>(d1.get()));
                }
    }
}


template<typename GridLayoutT>
void check_level(auto& rm, auto& level)
{
    using ResourcesManager_t = std::decay_t<decltype(rm)>;
    using FieldData_t        = typename ResourcesManager_t::UserField_t::patch_data_type;

    if constexpr (DEBUG_DIM == GridLayoutT::dimension)
    {
        auto& debugger = Debuggerino::INSTANCE();

        for (auto [name, info] : rm.all_resources())
            for (int i = 0; i < level.getLocalNumberOfPatches() - 1; ++i)
                for (int j = i + 1; j < level.getLocalNumberOfPatches(); ++j)
                {
                    auto const& p0 = *(level.getPatch(i));
                    auto const& p1 = *(level.getPatch(j));
                    auto const& d0 = p0.getPatchData(info.id);
                    auto const& d1 = p1.getPatchData(info.id);
                    auto f0        = dynamic_cast<FieldData_t*>(d0.get());
                    if (f0)
                        check_fields(f0, dynamic_cast<FieldData_t*>(d1.get()));
                }
    }
}

} // namespace PHARE::core

#if !PHARE_DEBUGGERINO

#define PHARE_DEBUG_SCOPE(key)
#define PHARE_DEBUG_SET(...)
#define PHARE_DEBUG_FIELDS(...)
#define PHARE_DEBUG_FIELD(...)
#define PHARE_DEBUG_ALL_FIELDS(...)
#define PHARE_DEBUG_CHECK_ALL_OVERLAPS(...)
#define PHARE_DEBUG_CHECK_LEVEL(...)

#endif // PHARE_DEBUGGERINO


#endif /*PHARE_CORE_DEBUG_HPP*/
