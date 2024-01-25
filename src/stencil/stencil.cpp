#include "stencil/stencil.hpp"
#include "slave/stencil_slave.hpp"
#include "stencil/bmp_image.hpp"
#include "stencil/boundary_matrix.hpp"

#include <algorithm>
#include <string>

void Stencil::initialize_matrix() {
    // We set off with an initial solution of 0. The left and right boundary are fixed at 1, while
    // the upper and lower boundaries are set to 0.
    matrix.fill_boundary(BoundaryPosition::Left, 1.f);
    matrix.fill_boundary(BoundaryPosition::Right, 1.f);

    // We alos initialize `result`, because we will swap `matrix` and `result` in the next
    // iteration and `result` will become the input.
    result.fill_boundary(BoundaryPosition::Left, 1.f);
    result.fill_boundary(BoundaryPosition::Right, 1.f);
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
    default:
        break;
    }
    auto const end = std::chrono::steady_clock::now();

    return end - start;
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
