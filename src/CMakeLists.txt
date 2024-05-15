add_subdirectory(stencil)

stencil_add_executable(
    stencil_main
    SOURCES main.cpp program_options.cpp
    LIBS stencil CLI11::CLI11
)
