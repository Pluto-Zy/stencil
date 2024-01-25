#ifndef STENCIL_SLAVE_STENCIL_SLAVE_HPP
#define STENCIL_SLAVE_STENCIL_SLAVE_HPP

#include <athread.h>
#include <crts.h>

#include <stencil/boundary_matrix.hpp>

#ifdef __sw_slave__
    #define SLAVE_FUN(name) name
#endif

struct Arguments {
    /// The side length of the block obtained by each CPE. If 64 CPEs are not enough to process the
    /// entire matrix in one round, then they divide the matrix into large chunks and process them
    /// in multiple rounds.
    unsigned block_size;
    /// Number of iterations.
    unsigned iterations;
    /// Points to the elements of the input matrix.
    BoundaryMatrixView<float> input;
    /// Points to the elements of the output matrix.
    BoundaryMatrixView<float> output;
};

/// Iteratively performs the stencil operation, where all CPEs use DMA to load elements from main
/// memory and additionally load required neighbor elements.
extern "C" void SLAVE_FUN(stencil_iterate_dma)(Arguments* args);
/// Iteratively performs the stencil operation, where all CPEs only load their own elements through
/// DMA, and use RMA to obtain neighbor elements from other CPEs.
extern "C" void SLAVE_FUN(stencil_iterate_rma)(Arguments* args);

#endif  // STENCIL_SLAVE_STENCIL_SLAVE_HPP
