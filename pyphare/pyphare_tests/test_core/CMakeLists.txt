cmake_minimum_required (VERSION 3.20.1)

project(test-pyphare-core)


add_python3_test(test-pyphare-box test_box.py ${PROJECT_SOURCE_DIR})
add_python3_test(test-pyphare-utilities test_phare_utilities.py ${PROJECT_SOURCE_DIR})
