add_subdirectory(slave)

stencil_add_library(
    stencil
    SOURCES stencil.cpp bmp_image.cpp program_options.cpp
    LIBS stencil_slave
)
