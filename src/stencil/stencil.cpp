#include "stencil/stencil.hpp"
#include "slave/stencil_slave.hpp"
#include "stencil/bmp_image.hpp"
#include "stencil/boundary_matrix.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <string_view>
#include <tuple>
#include <unordered_map>

void Stencil::initialize_matrix() {
    std::tie(matrix, result) = generate_initialized_matrix(
        options.matrix_size,
        options.matrix_size,
        options.radius,
        options.radius
    );
}

auto Stencil::run(Stencil::InputMethod method) -> std::chrono::steady_clock::duration {
    initialize_matrix();

    Arguments args = {
        /*block_size=*/options.block_size,
        /*iterations=*/options.iterations,
        /*input=*/matrix.borrow(),
        /*output=*/result.borrow(),
    };

    auto const start = std::chrono::steady_clock::now();
    switch (method) {
    case DMA:
        athread_spawn((void*) SLAVE_FUN(stencil_iterate_dma), &args);
        athread_join();
        break;
    case DMA_SLAVE_PACK:
        athread_spawn((void*) SLAVE_FUN(stencil_iterate_dma_slave_pack), &args);
        athread_join();
        break;
    case DMA_STATIC_UNROLL:
        athread_spawn((void*) SLAVE_FUN(stencil_iterate_dma_static_unroll), &args);
        athread_join();
        break;
    case RMA:
        athread_spawn((void*) SLAVE_FUN(stencil_iterate_rma), &args);
        athread_join();
        break;
    default:
        break;
    }
    auto const end = std::chrono::steady_clock::now();

    return end - start;
}

auto Stencil::run(std::string_view method_names) -> std::chrono::steady_clock::duration {
    // clang-format off
    static std::unordered_map<std::string_view, InputMethod> const method_map = {
        { "DMA", DMA },
        { "DMAStaticUnroll", DMA_STATIC_UNROLL },
        { "DMASlavePack", DMA_SLAVE_PACK },
        { "RMA", RMA },
        { "RMAWithoutSync", RMA },
    };
    // clang-format on

    auto const iter = method_map.find(method_names);
    assert(iter != method_map.end());

    return run(iter->second);
}

auto Stencil::check_result() const -> bool {
    // Generate input for the naive implementation.
    auto [naive_input, naive_output] = generate_initialized_matrix(
        options.matrix_size,
        options.matrix_size,
        options.radius,
        options.radius
    );

    // Run a naive implementation of the algorithm to produce correct results.
    float const avg = 1.f
        / static_cast<float>((naive_input.boundary_width() + naive_input.boundary_height()) * 2);

    // In the current implementation, the input buffer and output buffer are swapped every
    // iteration, which is a shallow swap. This makes it impossible to tell which buffer holds the
    // final result after the last iteration, so we use this variable to indicate the exchange
    // status.
    bool swapped = false;

    for (unsigned i = 0; i != options.iterations; ++i) {
        for (unsigned row = naive_input.boundary_height();
             row != naive_input.height_with_boundary() - naive_input.boundary_height();
             ++row)
        {
            for (unsigned col = naive_input.boundary_width();
                 col != naive_input.width_with_boundary() - naive_input.boundary_width();
                 ++col)
            {
                float sum = 0.f;

                // left
                for (unsigned nc = col - naive_input.boundary_width(); nc != col; ++nc) {
                    sum += naive_input.elem_with_boundary_at(row, nc);
                }

                // right
                for (unsigned nc = col; nc != col + naive_input.boundary_width(); ++nc) {
                    sum += naive_input.elem_with_boundary_at(row, nc + 1);
                }

                // top
                for (unsigned nr = row - naive_input.boundary_height(); nr != row; ++nr) {
                    sum += naive_input.elem_with_boundary_at(nr, col);
                }

                // bottom
                for (unsigned nr = row; nr != row + naive_input.boundary_height(); ++nr) {
                    sum += naive_input.elem_with_boundary_at(nr + 1, col);
                }

                naive_output.elem_with_boundary_at(row, col) = sum * avg;
            }
        }

        std::swap(naive_input, naive_output);
        swapped = !swapped;
    }

    // Check the result.
    BoundaryMatrix<float> const& compared = swapped ? result : matrix;
    for (unsigned row = 0; row != naive_input.height(); ++row) {
        for (unsigned col = 0; col != naive_input.width(); ++col) {
            if (std::fabs(naive_input.elem_at(row, col) - compared.elem_at(row, col)) >= 1e-4) {
                std::printf(
                    "invalid result at (%u, %u): %.15f vs %.15f\n",
                    row,
                    col,
                    naive_input.elem_at(row, col),
                    compared.elem_at(row, col)
                );
                return false;
            }
        }
    }

    return true;
}

auto Stencil::to_bmp() const -> BMPImage {
    std::vector<BMPImage::Pixel> image_data;
    image_data.reserve(options.matrix_size * options.matrix_size);

    for (unsigned i = 0; i != result.height(); ++i) {
        std::transform(
            &matrix.elem_at(i, 0),
            &matrix.elem_at(i, result.width()),
            std::back_inserter(image_data),
            [](float value) {
                BMPImage::Pixel rgb;
                if (value < 0.25) {
                    rgb.r = 0;
                    rgb.g = 4 * value * 255;
                    rgb.b = 255;
                } else if (value < 0.5) {
                    rgb.r = 0;
                    rgb.g = 255;
                    rgb.b = static_cast<unsigned char>((1 + 4 * (0.25 - value)) * 255);
                } else if (value < 0.75) {
                    rgb.r = static_cast<unsigned char>((4 * (value - 0.5)) * 255);
                    rgb.g = 255;
                    rgb.b = 0;
                } else {
                    assert(value <= 1.0);
                    rgb.r = 255;
                    rgb.g = static_cast<unsigned char>((1 + 4 * (0.75 - value)) * 255);
                    rgb.b = 0;
                }
                return rgb;
            }
        );
    }

    return BMPImage(std::move(image_data), options.matrix_size, options.matrix_size);
}

auto Stencil::generate_initialized_matrix(
    std::size_t width,
    std::size_t height,
    unsigned boundary_width,
    unsigned boundary_height
) -> std::tuple<BoundaryMatrix<float>, BoundaryMatrix<float>> {
    BoundaryMatrix<float> input1(width, height, boundary_width, boundary_height);
    BoundaryMatrix<float> input2(width, height, boundary_width, boundary_height);

    // We set off with an initial solution of 0. The left and right boundary are fixed at 1, while
    // the upper and lower boundaries are set to 0.
    input1.fill_boundary(BoundaryPosition::Left, 1.f);
    input1.fill_boundary(BoundaryPosition::Right, 1.f);
    input2.fill_boundary(BoundaryPosition::Left, 1.f);
    input2.fill_boundary(BoundaryPosition::Right, 1.f);

    return std::make_tuple(std::move(input1), std::move(input2));
}
