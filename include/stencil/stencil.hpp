#ifndef STENCIL_STENCIL_HPP
#define STENCIL_STENCIL_HPP

#include <chrono>
#include <string_view>

#include <stencil/bmp_image.hpp>
#include <stencil/boundary_matrix.hpp>
#include <stencil/program_options.hpp>

class Stencil {
public:
    /// Specifies the method for CPE to load input elements.
    enum InputMethod {
        /// All CPEs use DMA to load elements from main memory and additionally load required
        /// neighbor elements.
        DMA,
        /// All CPEs only load their own elements through DMA, and use RMA to obtain neighbor
        /// elements from other CPEs.
        RMA,
        /// Iteratively performs the stencil operation, where all CPEs use DMA to load elements
        /// from main memory and additionally load required neighbor elements. Local blocks of each
        /// CPE (including extra boundaries) are stored in contiguous memory. When exchanging
        /// non-contiguous memory areas, they will first pack them into contiguous memory areas.
        DMA_SLAVE_PACK,
        /// Iteratively performs the stencil operation, where all CPEs use DMA to load elements
        /// from main memory and additionally load required neighbor elements. The implementation
        /// is equivalent to the general implementation in `stencil_iterate_dma`, but enables loop
        /// unrolling at compile time by enumerating the neighborhood widths.
        DMA_STATIC_UNROLL,
    };

    Stencil() = default;
    explicit Stencil(ProgramOptions options) :
        options(options),
        matrix(options.matrix_size, options.matrix_size, options.radius, options.radius),
        result(options.matrix_size, options.matrix_size, options.radius, options.radius) { }

    void initialize_matrix();

    auto run(InputMethod method) -> std::chrono::steady_clock::duration;
    auto run(std::string_view method_name) -> std::chrono::steady_clock::duration;

    auto to_bmp() const -> BMPImage;

private:
    ProgramOptions options;
    /// The input of an iteration.
    BoundaryMatrix<float> matrix;
    /// The result of an iteration.
    BoundaryMatrix<float> result;
};

#endif  // STENCIL_STENCIL_HPP
