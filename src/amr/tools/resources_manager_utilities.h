#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_H
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_H
#include <type_traits>

namespace PHARE
{
template<typename...>
using tryToInstanciate = void;

template<typename ResourcesUser, typename Attempt = void>
struct has_field : std::false_type
{
};

template<typename ResourcesUser, typename Attempt = void>
struct has_particles : std::false_type
{
};

template<typename ResourcesUser, typename Attempt = void>
struct has_sub_resources : std::false_type
{
};

/** \brief has_field is a traits that permit to check if a ResourcesUser
 * has field
 *
 */
template<typename ResourcesUser>
struct has_field<ResourcesUser, tryToInstanciate<decltype(
                                    std::declval<ResourcesUser>().getFieldNamesAndQuantities())>>
    : std::true_type
{
};




/** \brief has_particles is a traits that permit to check if a ResourcesUser
 * has particles
 *
 */
template<typename ResourcesUser>
struct has_particles<ResourcesUser, tryToInstanciate<decltype(
                                        std::declval<ResourcesUser>().getParticleArrayNames())>>
    : std::true_type
{
};



/** \brief has_sub_resources is a traits that permit to check if a ResourcesUser
 * have other ResourcesUser
 *
 */
template<typename ResourcesUser>
struct has_sub_resources<ResourcesUser, tryToInstanciate<decltype(
                                            std::declval<ResourcesUser>().getSubResourcesObject())>>
    : std::true_type
{
};


struct dummy
{
    using type              = int;
    static const type value = 0;
};


/** \brief isFieldType will detect if ResourcesUser expose only Field */
template<typename ResourcesUser>
using isFieldType
    = std::enable_if_t<has_field<ResourcesUser>::value && !has_particles<ResourcesUser>::value
                           && !has_sub_resources<ResourcesUser>::value,
                       dummy::type>;



/** \brief isParticlesType will detect if ResourcesUser expose only particles */
template<typename ResourcesUser>
using isParticlesType
    = std::enable_if_t<!has_field<ResourcesUser>::value && has_particles<ResourcesUser>::value
                           && !has_sub_resources<ResourcesUser>::value,
                       dummy::type>;



/** \brief isSubResourcesType will detect if ResourcesUser expose only sub resources */
template<typename ResourcesUser>
using isSubResourcesType
    = std::enable_if_t<!has_field<ResourcesUser>::value && !has_particles<ResourcesUser>::value
                           && has_sub_resources<ResourcesUser>::value,
                       dummy::type>;



/** \brief isFieldAndParticlesType will detect if ResourcesUser expose only Field and Particles */
template<typename ResourcesUser>
using isFieldAndParticlesType
    = std::enable_if_t<has_field<ResourcesUser>::value && has_particles<ResourcesUser>::value
                           && !has_sub_resources<ResourcesUser>::value,
                       dummy::type>;



/** \brief isFieldAndSubResourcesType will detect if ResourcesUser expose only Field and sub
 * resources */
template<typename ResourcesUser>
using isFieldAndSubResourcesType
    = std::enable_if_t<has_field<ResourcesUser>::value && !has_particles<ResourcesUser>::value
                           && has_sub_resources<ResourcesUser>::value,
                       dummy::type>;



#if PARTICLES_AND_SUB_RESOURCES_EXIST
template<typename ResourcesUser>
using isParticlesAndSubResourcesType
    = std::enable_if_t<!has_field<ResourcesUser>::value && has_particles<ResourcesUser>::value
                           && has_sub_resources<ResourcesUser>::value,
                       std::true_type>;
#endif

#if FIELD_PARTICLES_AND_SUB_RESOURCES_EXIST
template<typename ResourcesUser>
using isAllType
    = std::enable_if_t<has_field<ResourcesUser>::value && has_particles<ResourcesUser>::value
                           && has_sub_resources<ResourcesUser>::value,
                       std::true_type>;
#endif

/** UseResourcePtr is used to select the resources patch data */
struct UseResourcePtr
{
};


/** UseNullPtr is used to select a nullptr with the correct type */
struct UseNullPtr
{
};


template<typename RequestedPtr>
using ifNullPtr = std::enable_if_t<
    std::is_same<typename std::remove_reference<RequestedPtr>::type, UseNullPtr>::value,
    dummy::type>;


template<typename RequestedPtr>
using ifResourcePtr = std::enable_if_t<
    std::is_same<typename std::remove_reference<RequestedPtr>::type, UseResourcePtr>::value,
    dummy::type>;



} // namespace PHARE

#endif // PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_H
