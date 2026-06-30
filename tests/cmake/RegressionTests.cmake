function(add_regression_test TEST_NAME CONFIG_FILE OUTPUT_FILE)
    set(CONFIG_PATH "${TEST_CONFIG_DIR}/${CONFIG_FILE}")
    set(REFERENCE_PATH "${TEST_REFERENCE_DIR}/${OUTPUT_FILE}")
    set(OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_FILE}")

    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND} -E env
            ${CMAKE_COMMAND}
            -DMONOFONIC_EXECUTABLE=$<TARGET_FILE:monofonIC>
            -DCONFIG_FILE=${CONFIG_PATH}
            -DOUTPUT_FILE=${OUTPUT_PATH}
            -DREFERENCE_FILE=${REFERENCE_PATH}
            -DCOMPARE_SCRIPT=${COMPARE_SCRIPT}
            -DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT 300
        LABELS "regression"
    )

    if(NOT EXISTS ${REFERENCE_PATH})
        set_tests_properties(${TEST_NAME} PROPERTIES DISABLED TRUE)
        message(STATUS "Test ${TEST_NAME}: Reference file not found - test will be disabled")
        message(STATUS "  Generate references with: cd build/tests && bash ${TEST_SCRIPT_DIR}/generate_references.sh")
    else()
        message(STATUS "Test ${TEST_NAME}: Registered")
    endif()
endfunction()

add_regression_test(test_1lpt_sc_generic test_1lpt_sc_generic.conf test_1lpt_sc_generic.hdf5)
add_regression_test(test_2lpt_sc_gadget test_2lpt_sc_gadget.conf test_2lpt_sc_gadget.hdf5)
add_regression_test(test_3lpt_bcc_swift test_3lpt_bcc_swift.conf test_3lpt_bcc_swift.hdf5)
add_regression_test(test_2lpt_baryons_generic test_2lpt_baryons_generic.conf test_2lpt_baryons_generic.hdf5)
add_regression_test(test_2lpt_baryons_vrel_gadget test_2lpt_baryons_vrel_gadget.conf test_2lpt_baryons_vrel_gadget.hdf5)

add_test(NAME compare_script_help COMMAND ${Python3_EXECUTABLE} ${COMPARE_SCRIPT} --help)
set_tests_properties(compare_script_help PROPERTIES TIMEOUT 10 LABELS "sanity")
