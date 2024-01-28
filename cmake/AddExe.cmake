function(stencil_add_library name)
    cmake_parse_arguments(LIB "" "FOLDER" "SOURCES;LIBS" ${ARGN})
    add_library(${name} OBJECT ${LIB_SOURCES})

    target_compile_options(${name} PRIVATE
            -Wall -Wextra -Wpedantic -Wno-unused-function
            -mhost -mieee)

    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${name} PUBLIC ${LIB_LIBS})

    foreach(library ${LIB_LIBS})
        get_target_property(library_type ${library} TYPE)
        if (library_type STREQUAL "OBJECT_LIBRARY")
            target_link_libraries(${name} PUBLIC $<TARGET_OBJECTS:${library}>)
        endif()
    endforeach()

    if (LIB_FOLDER)
        set_target_properties(${name} PROPERTIES FOLDER ${LIB_FOLDER})
    endif ()
endfunction(stencil_add_library)

function(stencil_add_executable name)
    cmake_parse_arguments(EXE "" "FOLDER" "SOURCES;LIBS" ${ARGN})
    add_executable(${name} ${EXE_SOURCES})

    target_compile_options(${name} PRIVATE
            -Wall -Wextra -Wpedantic -Wno-unused-function
            -mhost -mieee)

    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${name} PUBLIC ${EXE_LIBS})

    # foreach(library ${EXE_LIBS})
    #     get_target_property(library_type ${library} TYPE)
    #     if (library_type STREQUAL "OBJECT_LIBRARY")
    #         target_link_libraries(${name} PUBLIC $<TARGET_OBJECTS:${library}>)
    #     endif()
    # endforeach()

    if (EXE_FOLDER)
        set_target_properties(${name} PROPERTIES FOLDER ${EXE_FOLDER})
    endif ()
endfunction(stencil_add_executable)

function(stencil_add_slave_library name)
    cmake_parse_arguments(LIB "" "" "SOURCES;LIBS" ${ARGN})
    add_library(${name} OBJECT ${LIB_SOURCES})

    target_compile_options(${name} PRIVATE
            -Wall -Wextra -Wpedantic
            -mslave -mieee -msimd
            -funroll-all-loops)

    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include/experiment/slave)
    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${name} PUBLIC ${LIB_LIBS})
    target_link_options(${name} PUBLIC -mhybrid)

    # foreach(library ${LIB_LIBS})
    #     get_target_property(library_type ${library} TYPE)
    #     if (library_type STREQUAL "OBJECT_LIBRARY")
    #         target_link_libraries(${name} PUBLIC $<TARGET_OBJECTS:${library}>)
    #     endif()
    # endforeach()
endfunction(stencil_add_slave_library)
