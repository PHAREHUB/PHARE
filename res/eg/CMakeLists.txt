
# considered heavy, so as 3d
if(test AND testMPI AND executeNotebooks AND HighFive AND PHARE_PYTEST_SIMULATORS GREATER_EQUAL 3)

    function(convert_notebook_to_script inpath outpath)
      execute_process(
        COMMAND ./convert_to_script.sh "${inpath}" "${outpath}"
        COMMAND_ERROR_IS_FATAL ANY
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    endfunction(convert_notebook_to_script)

    convert_notebook_to_script("1d/jupyter/strong/strong.ipynb" "${CMAKE_CURRENT_BINARY_DIR}/1d/jupyter/strong/strong_notebook")
    convert_notebook_to_script("1d/jupyter/weak/weak.ipynb" "${CMAKE_CURRENT_BINARY_DIR}/1d/jupyter/weak/weak_notebook")

    # internally these tests (can) call mpirun - so we pretend they are serial here
    phare_mpi_python3_exec(9 1 notebook_strong strong_notebook.py "${CMAKE_CURRENT_BINARY_DIR}/1d/jupyter/strong")
    phare_mpi_python3_exec(9 1 notebook_weak weak_notebook.py "${CMAKE_CURRENT_BINARY_DIR}/1d/jupyter/weak")

endif(test AND testMPI AND executeNotebooks AND HighFive AND PHARE_PYTEST_SIMULATORS GREATER_EQUAL 3)
