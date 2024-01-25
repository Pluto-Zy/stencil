#include <climits>
#include <cstring>
#include <fstream>
#include <stencil/bmp_image.hpp>
#include <string_view>

void BMPImage::save(std::string_view filename) const {
    std::ofstream out(filename.data(), std::ios::binary);
    assert(out);

    std::uint32_t const padding_size = (4 - (width * sizeof(Pixel)) % 4) % 4;

    // Write file header.
    write_file_header(out, padding_size);
    // Write info header.
    write_info_header(out);

    // Write image data.
    unsigned char const padding[4] {};
    for (std::uint32_t i = 0; i != height; ++i) {
        out.write(reinterpret_cast<char const*>(&image_data[i * width]), width * sizeof(Pixel));

        // Write padding.
        out.write(reinterpret_cast<char const*>(padding), padding_size);
    }
}

void BMPImage::write_file_header(std::ofstream& out, std::uint32_t padding_size) const {
    unsigned char header[] = {
        'B', 'M',  // signature
        0,   0,   0, 0,  // image file size in bytes
        0,   0,   0, 0,  // reserved
        54,  0,   0, 0,  // start of pixel array
    };

    std::uint32_t const width_in_bytes = width * sizeof(Pixel) + padding_size;
    std::uint32_t const file_size = 54 + width_in_bytes * height;

    // Write file size to the header.
    std::memcpy(header + 2, &file_size, sizeof(std::uint32_t));

    // Write the header to file.
    out.write(reinterpret_cast<char*>(header), sizeof(header));
}

void BMPImage::write_info_header(std::ofstream& out) const {
    // clang-format off
    unsigned char header[40] = {
        40, 0, 0, 0,  // header size
        0, 0, 0, 0,  // image width
        0, 0, 0, 0,  // image height
        1, 0,  // number of color planes
        sizeof(Pixel) * CHAR_BIT, 0,  // bits per pixel
    };
    // clang-format on

    // Write image size to the header.
    std::uint32_t width = this->width, height = this->height;
    std::memcpy(header + 4, &width, sizeof(width));
    std::memcpy(header + 8, &height, sizeof(height));

    // Write the header to file.
    out.write(reinterpret_cast<char*>(header), sizeof(header));
}
