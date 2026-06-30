if(MPI_CXX_FOUND)
    set(MPI_TEST_SCRIPT "${TEST_SCRIPT_DIR}/test_mpi_consistency.sh")
    set(MPI_TEST_CONFIG "${TEST_CONFIG_DIR}/test_2lpt_mpi.conf")
    add_test(
        NAME test_mpi_consistency
        COMMAND bash ${MPI_TEST_SCRIPT}
            $<TARGET_FILE:monofonIC>
            ${MPI_TEST_CONFIG}
            ${COMPARE_SCRIPT}
            ${Python3_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(test_mpi_consistency PROPERTIES
        TIMEOUT 600
        LABELS "mpi;regression"
    )
    message(STATUS "Test test_mpi_consistency: Registered (tests with 1, 2, and 4 MPI tasks)")
else()
    message(STATUS "Test test_mpi_consistency: Skipped (MPI not available)")
endif()
