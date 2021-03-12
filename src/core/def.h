
#ifndef PHARE_CORE_DEF_H
#define PHARE_CORE_DEF_H


#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_DEBUG_DO(...) __VA_ARGS__

template<class... Args>
void PHARE_DEBUG_PRINT(Args... args)
{
    (std::cout << ... << args) << "\n";
}
#else
#define PHARE_DEBUG_DO(...)

template<class... Args>
void PHARE_DEBUG_PRINT(Args... args)
{
}
#endif

#endif // PHARE_CORE_DEF_H
