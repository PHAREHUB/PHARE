#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#define NO_DISCARD [[nodiscard]]


#define _PHARE_TO_STR(x) #x // convert macro text to string
#define PHARE_TO_STR(x) _PHARE_TO_STR(x)

#define PHARE_TOKEN_PASTE(x, y) x##y
#define PHARE_STR_CAT(x, y) PHARE_TOKEN_PASTE(x, y)

#endif // PHARE_CORE_DEF_HPP
