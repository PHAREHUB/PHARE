#ifndef PHARE_BENCH_CORE_NO_VECTORIZE_HPP
#define PHARE_BENCH_CORE_NO_VECTORIZE_HPP

// keeps a function scalar on gcc, for baselines: "#pragma GCC novector" is gcc >= 14 only
//  clang has no per function equivalent, use "#pragma clang loop vectorize(disable)"
#if !defined(__clang__) && defined(__GNUC__)
#define PHARE_BENCH_NO_VECTORIZE __attribute__((optimize("no-tree-vectorize")))
#else
#define PHARE_BENCH_NO_VECTORIZE
#endif

#endif /* PHARE_BENCH_CORE_NO_VECTORIZE_HPP */
