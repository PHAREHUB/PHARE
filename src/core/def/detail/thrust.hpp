#ifndef PHARE_CORE_DEF_THRUST_HPP
#define PHARE_CORE_DEF_THRUST_HPP




// clang-format off
#if PHARE_HAVE_GPU && !defined(PHARE_HAVE_THRUST) && __has_include(<thrust/iterator/zip_iterator.h>)
  #define PHARE_HAVE_THRUST 1
#endif // PHARE_HAVE_THRUST

#if !defined(PHARE_HAVE_THRUST)
  #define PHARE_HAVE_THRUST 0
#endif


#if PHARE_HAVE_THRUST

  #include <thrust/sort.h>
  #include <thrust/execution_policy.h>
  #include <thrust/iterator/zip_iterator.h>

  #define PHARE_WITH_THRUST(...) __VA_ARGS__
  #define PHARE_WITH_THRUST_ELSE(...)
  #define PHARE_WITH_THRUST_ELSE_THROW(...) __VA_ARGS__

#else // !__has_include(<thrust/iterator/zip_iterator.h>)

  #define PHARE_HAVE_THRUST 0
  #define PHARE_WITH_THRUST(...)
  #define PHARE_WITH_THRUST_ELSE(...) __VA_ARGS__
  #define PHARE_WITH_THRUST_ELSE_THROW(...) throw std::runtime_error("Thrust not found!");

#endif // PHARE_HAVE_THRUST

// clang-format on

#endif /* PHARE_CORE_DEF_THRUST_HPP */
