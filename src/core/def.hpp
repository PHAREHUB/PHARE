#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP



#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_DEBUG_DO(...) __VA_ARGS__
#else
#define PHARE_DEBUG_DO(...)
#endif


#define _PHARE_TO_STR(x) #x // convert macro text to string
#define PHARE_TO_STR(x) _PHARE_TO_STR(x)



#endif /* PHARE_CORE_DEF_HPP */
