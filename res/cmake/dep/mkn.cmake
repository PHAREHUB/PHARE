# !!!!!!

if(withMkn)

    phare_github_get_or_update(mkn.kul ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mkn/kul mkn/mkn.kul master)
    phare_github_get_or_update(mkn.gpu ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mkn/gpu mkn/mkn.gpu master)

    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mkn/kul/inc)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/subprojects/mkn/gpu/inc)

    # phare_github_get_or_update(cccl ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/nvcccl NVIDIA/cccl main)
    # add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/subprojects/nvcccl)
    # set (PHARE_BASE_LIBS ${PHARE_BASE_LIBS} CUB::CUB Thrust::Thrust)

endif(withMkn)
