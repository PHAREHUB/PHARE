#ifndef PHARE_MAGNETIC_FIELD_COARSENER
#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
     * given coarse node
     *
     * A MagneticFieldCoarsener object is created each time the refine() method of the
     * FieldCoarsenOperator is called and its operator() is called for each coarse index.
     * It is the default coarsening policy and used for any field that does not come with
     * specific constraints (such as conserving some property in the coarsening process).
     */
    template<std::size_t dimension>
    class MagneticFieldCoarsener
    {
    public:
        MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                               SAMRAI::hier::Box const& sourceBox,
                               SAMRAI::hier::Box const& destinationBox,
                               SAMRAI::hier::IntVector const& ratio)
        {
        }
        template<typename FieldT>
        void operator()(FieldT const& fineField, FieldT& coarseField,
                        core::Point<int, dimension> coarseIndex)
        {
        }
    };
} // namespace amr
} // namespace PHARE




#ifndef PHARE_MAGNETIC_FIELD_COARSENER
#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
     * given coarse node
     *
     * A MagneticFieldCoarsener object is created each time the refine() method of the
     * FieldCoarsenOperator is called and its operator() is called for each coarse index.
     * It is the default coarsening policy and used for any field that does not come with
     * specific constraints (such as conserving some property in the coarsening process).
     */
    template<std::size_t dimension>
    class MagneticFieldCoarsener
    {
    public:
        MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                               SAMRAI::hier::Box const& sourceBox,
                               SAMRAI::hier::Box const& destinationBox,
                               SAMRAI::hier::IntVector const& ratio)
        {
        }
        template<typename FieldT>
        void operator()(FieldT const& fineField, FieldT& coarseField,
                        core::Point<int, dimension> coarseIndex)
        {
            // For the moment we only take the case of field with the same centering
            TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

            core::Point<int, dimension> fineStartIndex;

            fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
            coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


            double coarseValue = 0.;
        }
    };
} // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
     * given coarse node
     *
     * A MagneticFieldCoarsener object is created each time the refine() method of the
     * FieldCoarsenOperator is called and its operator() is called for each coarse index.
     * It is the default coarsening policy and used for any field that does not come with
     * specific constraints (such as conserving some property in the coarsening process).
     */
    template<std::size_t dimension>
    class MagneticFieldCoarsener
    {
    public:
        MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                               SAMRAI::hier::Box const& sourceBox,
                               SAMRAI::hier::Box const& destinationBox,
                               SAMRAI::hier::IntVector const& ratio)
        {
        }
        template<typename FieldT>
        void operator()(FieldT const& fineField, FieldT& coarseField,
                        core::Point<int, dimension> coarseIndex)
        {
            // For the moment we only take the case of field with the same centering
            TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

            core::Point<int, dimension> fineStartIndex;
            if consteconstexpr

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
            coarseIndex = AMRToLocal(coarseIndex, destinationBox_);


            double coarseValue = 0.;
        }
    };
} // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
     * given coarse node
     *
     * A MagneticFieldCoarsener object is created each time the refine() method of the
     * FieldCoarsenOperator is called and its operator() is called for each coarse index.
     * It is the default coarsening policy and used for any field that does not come with
     * specific constraints (such as conserving some property in the coarsening process).
     */
    template<std::size_t dimension>
    class MagneticFieldCoarsener
    {
    public:
        MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                               SAMRAI::hier::Box const& sourceBox,
                               SAMRAI::hier::Box const& destinationBox,
                               SAMRAI::hier::IntVector const& ratio)
        {
        }
        template<typename FieldT>
        void operator()(FieldT const& fineField, FieldT& coarseField,
                        core::Point<int, dimension> coarseIndex)
        {
            // For the moment we only take the case of field with the same centering
            TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

            core::Point<int, dimension> fineStartIndex;
            if contexpr (dimension == 1)

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
            coarseIndex = AMRToLocal(coarseIndex, destinationBox_);


            double coarseValue = 0.;
        }
    };
} // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
     * given coarse node
     *
     * A MagneticFieldCoarsener object is created each time the refine() method of the
     * FieldCoarsenOperator is called and its operator() is called for each coarse index.
     * It is the default coarsening policy and used for any field that does not come with
     * specific constraints (such as conserving some property in the coarsening process).
     */
    template<std::size_t dimension>
    class MagneticFieldCoarsener
    {
    public:
        MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                               SAMRAI::hier::Box const& sourceBox,
                               SAMRAI::hier::Box const& destinationBox,
                               SAMRAI::hier::IntVector const& ratio)
        {
        }
        template<typename FieldT>
        void operator()(FieldT const& fineField, FieldT& coarseField,
                        core::Point<int, dimension> coarseIndex)
        {
            // For the moment we only take the case of field with the same centering
            TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

            core::Point<int, dimension> fineStartIndex;
            if contexpr (dimension == 1)
            {
                ]

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;
            }
        };
    } // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio)
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;
            }
        };
    } // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio)
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;
            }
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio)
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;
            }
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio)
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ];
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;
            }
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY]) = 0.5*(fineField(fineStartIndex[dirX], fineStartIndex[dirY)
                                                                     + 2);
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is -dual-, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                // Bx is primal-dual, take average in Y
                if (true)
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                // Bx is primal-dual, take average in Y
                if (coarseField)
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                // Bx is primal-dual, take average in Y
                if (coarseField.physicalQuantity() == hybricore::HybridQuantity)
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if (centering)
                    // Bx is primal-dual, take average in Y
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if constexpr (dimension = +1)
                    // Bx is primal-dual, take average in Y
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if constexpr (dimension == 2)
                {
                    // Bx is primal-dual, take average in Y
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                    // By is dual-primal, take average in X
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                    // Bz is dual-dual, take average in XY
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.25
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
                }
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if constexpr (dimension == 2)
                {
                    if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                                  and centering[idirY] == core::QtyCentering::dual)
                        // Bx is primal-dual, take average in Y
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));

                    // By is dual-primal, take average in X
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                    // Bz is dual-dual, take average in XY
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.25
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
                }
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if constexpr (dimension == 2)
                {
                    if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                                  and centering[idirY] == core::QtyCentering::dual)
                    {
                        // Bx is primal-dual, take average in Y
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
                    }
                    // By is dual-primal, take average in X
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));


                    // Bz is dual-dual, take average in XY
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.25
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
                }
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if constexpr (dimension == 2)
                {
                    if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                                  and centering[idirY] == core::QtyCentering::dual)
                    {
                        // Bx is primal-dual, take average in Y
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
                    }
                    if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                                  and centering[idirY] == core::QtyCentering::dual)
                    {
                        // By is dual-primal, take average in X
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
                    }
                    // Bz is dual-dual, take average in XY
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.25
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
                }
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




namespace PHARE
{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if constexpr (dimension == 2)
                {
                    if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                                  and centering[idirY] == core::QtyCentering::dual)
                    {
                        // Bx is primal-dual, take average in Y
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
                    }
                    if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                                  and centering[idirY] == core::QtyCentering::dual)
                    {
                        // By is dual-primal, take average in X
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
                    }
                    // Bz is dual-dual, take average in XY
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.25
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
                }
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




{namespace amr{using core::dirX;
using core::dirY;
using core::dirZ;
/** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
 * a given coarse node
 *
 * A MagneticFieldCoarsener object is created each time the refine() method of the
 * FieldCoarsenOperator is called and its operator() is called for each coarse index.
 * It is the default coarsening policy and used for any field that does not come with
 * specific constraints (such as conserving some property in the coarsening process).
 */
template<std::size_t dimension>
class MagneticFieldCoarsener
{
public:
    MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                           SAMRAI::hier::Box const& sourceBox,
                           SAMRAI::hier::Box const& destinationBox,
                           SAMRAI::hier::IntVector const& ratio),
        : sourceBox_{sourceBox}, destinationBox_{destinationBox}
    {
    }
    template<typename FieldT>
    void operator()(FieldT const& fineField, FieldT& coarseField,
                    core::Point<int, dimension> coarseIndex)
    {
        // For the moment we only take the case of field with the same centering
        TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

        core::Point<int, dimension> fineStartIndex;
        if contexpr (dimension == 1)
        {
            fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
        }
        else if constexpr (dimension > 1)
        {
            fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
            if constepr (dimension > 2)
            {
                fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
            }
        }

        fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
        coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


        double coarseValue = 0.;

        // we're not supposed to know B has this specific centering here
        // hard coded for now but should instead used the layout to ask for centering

        if constexpr (dimension == 2)
        {
            if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                          and centering[idirY] == core::QtyCentering::dual)
            {
                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                          and centering[idirY] == core::QtyCentering::primal)
            {
                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            // Bz is dual-dual, take average in XY
            coarseField(coarseIndex[dirX], coarseIndex[dirY])
                = 0.25
                  * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                     + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                     + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                     + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
        }
    }
    SAMRAI::hier::Box const sourceBox_;
    SAMRAI::hier::Box const destinationBox_;
    int constexpr ratio_ = 2;
};
} // namespace amr
} // namespace PHARE




{
    namespace amr
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;
        /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
         * a given coarse node
         *
         * A MagneticFieldCoarsener object is created each time the refine() method of the
         * FieldCoarsenOperator is called and its operator() is called for each coarse index.
         * It is the default coarsening policy and used for any field that does not come with
         * specific constraints (such as conserving some property in the coarsening process).
         */
        template<std::size_t dimension>
        class MagneticFieldCoarsener
        {
        public:
            MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                                   SAMRAI::hier::Box const& sourceBox,
                                   SAMRAI::hier::Box const& destinationBox,
                                   SAMRAI::hier::IntVector const& ratio),
                : sourceBox_{sourceBox}, destinationBox_{destinationBox}
            {
            }
            template<typename FieldT>
            void operator()(FieldT const& fineField, FieldT& coarseField,
                            core::Point<int, dimension> coarseIndex)
            {
                // For the moment we only take the case of field with the same centering
                TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

                core::Point<int, dimension> fineStartIndex;
                if contexpr (dimension == 1)
                {
                    fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
                }
                else if constexpr (dimension > 1)
                {
                    fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                    if constepr (dimension > 2)
                    {
                        fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                    }
                }

                fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
                coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


                double coarseValue = 0.;

                // we're not supposed to know B has this specific centering here
                // hard coded for now but should instead used the layout to ask for centering

                if constexpr (dimension == 2)
                {
                    if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                                  and centering[idirY] == core::QtyCentering::dual)
                    {
                        // Bx is primal-dual, take average in Y
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
                    }
                    if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                                  and centering[idirY] == core::QtyCentering::primal)
                    {
                        // By is dual-primal, take average in X
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.5
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
                    }
                    if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                                  and centering[idirY] == core::QtyCentering::primal)
                    {
                        // Bz is dual-dual, take average in XY
                        coarseField(coarseIndex[dirX], coarseIndex[dirY])
                            = 0.25
                              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
                    }
                }
            }
            SAMRAI::hier::Box const sourceBox_;
            SAMRAI::hier::Box const destinationBox_;
            int constexpr ratio_ = 2;
        };
    } // namespace amr
} // namespace PHARE




using core::dirX;
using core::dirY;
using core::dirZ;
/** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
 * a given coarse node
 *
 * A MagneticFieldCoarsener object is created each time the refine() method of the
 * FieldCoarsenOperator is called and its operator() is called for each coarse index.
 * It is the default coarsening policy and used for any field that does not come with
 * specific constraints (such as conserving some property in the coarsening process).
 */
template<std::size_t dimension>
class MagneticFieldCoarsener
{
public:
    MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                           SAMRAI::hier::Box const& sourceBox,
                           SAMRAI::hier::Box const& destinationBox,
                           SAMRAI::hier::IntVector const& ratio),
        : sourceBox_{sourceBox}, destinationBox_{destinationBox}
    {
    }
    template<typename FieldT>
    void operator()(FieldT const& fineField, FieldT& coarseField,
                    core::Point<int, dimension> coarseIndex)
    {
        // For the moment we only take the case of field with the same centering
        TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

        core::Point<int, dimension> fineStartIndex;
        if contexpr (dimension == 1)
        {
            fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
        }
        else if constexpr (dimension > 1)
        {
            fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
            if constepr (dimension > 2)
            {
                fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
            }
        }

        fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
        coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


        double coarseValue = 0.;

        // we're not supposed to know B has this specific centering here
        // hard coded for now but should instead used the layout to ask for centering

        if constexpr (dimension == 2)
        {
            if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                          and centering[idirY] == core::QtyCentering::dual)
            {
                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                          and centering[idirY] == core::QtyCentering::primal)
            {
                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                          and centering[idirY] == core::QtyCentering::dual)
            {
                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
        }
    }
    SAMRAI::hier::Box const sourceBox_;
    SAMRAI::hier::Box const destinationBox_;
    int constexpr ratio_ = 2;
};
} // namespace amr
} // namespace PHARE




using core::dirX;
using core::dirY;
using core::dirZ;
/** @brief This class gives an operator() that performs the coarsening of N fine nodes onto
 * a given coarse node
 *
 * A MagneticFieldCoarsener object is created each time the refine() method of the
 * FieldCoarsenOperator is called and its operator() is called for each coarse index.
 * It is the default coarsening policy and used for any field that does not come with
 * specific constraints (such as conserving some property in the coarsening process).
 */
template<std::size_t dimension>
class MagneticFieldCoarsener
{
public:
    MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                           SAMRAI::hier::Box const& sourceBox,
                           SAMRAI::hier::Box const& destinationBox,
                           SAMRAI::hier::IntVector const& ratio),
        : sourceBox_{sourceBox}, destinationBox_{destinationBox}
    {
    }
    template<typename FieldT>
    void operator()(FieldT const& fineField, FieldT& coarseField,
                    core::Point<int, dimension> coarseIndex)
    {
        // For the moment we only take the case of field with the same centering
        TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

        core::Point<int, dimension> fineStartIndex;
        if contexpr (dimension == 1)
        {
            fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
        }
        else if constexpr (dimension > 1)
        {
            fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
            if constepr (dimension > 2)
            {
                fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
            }
        }

        fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
        coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


        double coarseValue = 0.;

        // we're not supposed to know B has this specific centering here
        // hard coded for now but should instead used the layout to ask for centering

        if constexpr (dimension == 2)
        {
            if constexpr (centering[idirX] == core::QtyQtyCentering::primal
                          and centering[idirY] == core::QtyCentering::dual)
            {
                // Bx is primal-dual, take average in Y
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                          and centering[idirY] == core::QtyCentering::primal)
            {
                // By is dual-primal, take average in X
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            if constexpr (centering[idirX] == core::QtyQtyCentering::dual
                          and centering[idirY] == core::QtyCentering::dual)
            {
                // Bz is dual-dual, take average in XY
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
            }
        }
    }
    SAMRAI::hier::Box const sourceBox_;
    SAMRAI::hier::Box const destinationBox_;
    int constexpr ratio_ = 2;
};
} // namespace amr
} // namespace PHARE




#endif // !PHARE_MAGNETIC_FIELD_COARSENER
