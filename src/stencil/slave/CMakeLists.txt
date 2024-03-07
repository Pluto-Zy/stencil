stencil_add_slave_library(
    stencil_slave
    SOURCES
    stencil_dma.cpp stencil_dma_static_unroll.cpp
    stencil_dma_slave_pack.cpp stencil_rma.cpp
)
