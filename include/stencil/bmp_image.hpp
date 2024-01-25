#ifndef STENCIL_BMP_IMAGE_HPP
#define STENCIL_BMP_IMAGE_HPP

#include <cassert>
#include <cstdint>
#include <fstream>
#include <string_view>
#include <vector>

class BMPImage {
public:
    struct Pixel {
        std::uint8_t b, g, r;
    };

    BMPImage() : width(0), height(0) { }
    explicit BMPImage(std::vector<Pixel> data, std::uint32_t width, std::uint32_t height) :
        image_data(std::move(data)), width(width), height(height) {
        assert(width * height == image_data.size());
    }

    void save(std::string_view filename) const;

private:
    std::vector<Pixel> image_data;
    std::uint32_t width;
    std::uint32_t height;

    void write_file_header(std::ofstream& out, std::uint32_t padding_size) const;
    void write_info_header(std::ofstream& out) const;
};

#endif  // STENCIL_BMP_IMAGE_HPP
