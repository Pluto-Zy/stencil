add_subdirectory(stencil)

stencil_add_executable(
    stencil_main
    SOURCES main.cpp
    LIBS stencil
)
