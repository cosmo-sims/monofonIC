function(add_plane_wave_test TEST_NAME CONFIG_FILE DIAG_FILENAME OUTPUT_FILENAME CHECK_SCRIPT)
    set(CONFIG_SOURCE "${TEST_CONFIG_DIR}/${CONFIG_FILE}")
    set(CONFIG_BINARY "${CMAKE_CURRENT_BINARY_DIR}/${CONFIG_FILE}")
    configure_file(${CONFIG_SOURCE} ${CONFIG_BINARY} COPYONLY)

    get_filename_component(CONFIG_DIR "${CONFIG_BINARY}" DIRECTORY)
    get_filename_component(CONFIG_NAME_WE "${CONFIG_BINARY}" NAME_WE)
    set(CONFIG_BASENAME "${CONFIG_DIR}/${CONFIG_NAME_WE}")
    set(DIAG_PATH "${CONFIG_BASENAME}_${DIAG_FILENAME}")
    set(OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_FILENAME}")
    set(CHECK_SCRIPT_PATH "${TEST_SCRIPT_DIR}/${CHECK_SCRIPT}")

    add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_COMMAND}
            -DMONOFONIC_EXECUTABLE=$<TARGET_FILE:monofonIC>
            -DCONFIG_FILE=${CONFIG_BINARY}
            -DOUTPUT_FILE=${OUTPUT_PATH}
            -DDIAG_FILE=${DIAG_PATH}
            -DCHECK_SCRIPT=${CHECK_SCRIPT_PATH}
            -DPYTHON_EXECUTABLE=${Python3_EXECUTABLE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_plane_wave_test.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT 300
        LABELS "regression;plane-wave"
    )
    message(STATUS "Test ${TEST_NAME}: Registered")
endfunction()

add_plane_wave_test(test_plane_wave_model1 test_plane_wave_model1.conf
    plane_wave_model1_diag.hdf5 plane_wave_model1_output.hdf5 test_plane_wave_model1.py)
add_plane_wave_test(test_plane_wave_model2 test_plane_wave_model2.conf
    plane_wave_model2_diag.hdf5 plane_wave_model2_output.hdf5 test_plane_wave_model2.py)
add_plane_wave_test(test_plane_wave_model3 test_plane_wave_model3.conf
    plane_wave_model3_diag.hdf5 plane_wave_model3_output.hdf5 test_plane_wave_model3.py)
