#ifndef STENCIL_STENCIL_HPP
#define STENCIL_STENCIL_HPP

#include <chrono>

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
    };

    Stencil() = default;
    explicit Stencil(ProgramOptions options) :
        options(options),
        matrix(
            options.matrix_size,
            options.matrix_size,
            options.neighbor_width,
            options.neighbor_width
        ),
        result(
            options.matrix_size,
            options.matrix_size,
            options.neighbor_width,
            options.neighbor_width
        ) { }

    void initialize_matrix();

    auto run(InputMethod method) -> std::chrono::steady_clock::duration;

    auto to_bmp() const -> BMPImage;

private:
    ProgramOptions options;
    /// The input of an iteration.
    BoundaryMatrix<float> matrix;
    /// The result of an iteration.
    BoundaryMatrix<float> result;
};

#endif  // STENCIL_STENCIL_HPP
