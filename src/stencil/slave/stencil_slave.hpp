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

extern "C" {
/// Iteratively performs the stencil operation, where all CPEs use DMA to load elements from main
/// memory and additionally load required neighbor elements.
void SLAVE_FUN(stencil_iterate_dma)(Arguments* args);
/// Iteratively performs the stencil operation, where all CPEs use DMA to load elements from main
/// memory and additionally load required neighbor elements. The implementation is equivalent to
/// the general implementation in `stencil_iterate_dma`, but enables loop unrolling at compile time
/// by enumerating the neighborhood widths.
void SLAVE_FUN(stencil_iterate_dma_static_unroll)(Arguments* args);
/// Iteratively performs the stencil operation, where all CPEs use DMA to load elements from main
/// memory and additionally load required neighbor elements. Local blocks of each CPE (including
/// extra boundaries) are stored in contiguous memory. When exchanging non-contiguous memory areas,
/// they will first pack them into contiguous memory areas.
///
/// This implementation is used to test the performance of packed and unpacked (compared to
/// stencil_iterate_dma).
void SLAVE_FUN(stencil_iterate_dma_slave_pack)(Arguments* args);
/// Iteratively performs the stencil operation, where all CPEs only load their own elements through
/// DMA, and use RMA to obtain neighbor elements from other CPEs.
void SLAVE_FUN(stencil_iterate_rma)(Arguments* args);
}

#endif  // STENCIL_SLAVE_STENCIL_SLAVE_HPP
