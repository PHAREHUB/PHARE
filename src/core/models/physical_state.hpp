#ifndef PHARE_PHYSICAL_STATE_H
#define PHARE_PHYSICAL_STATE_H




/**
 * @brief The PhysicalState class is an interface for concrete states holding data associated with a
 * concrete IPhysicalModel
 */
namespace PHARE
{
namespace core
{
    enum class Model { MHD, Hybrid };

    class PhysicalStateInitializer
    {
    };


    class IPhysicalState
    {
    };

} // namespace core
} // namespace PHARE



#endif // PHYSICAL_STATE_H
