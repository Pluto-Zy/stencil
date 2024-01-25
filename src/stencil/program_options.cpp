#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <unistd.h>

#include <stencil/program_options.hpp>

auto ProgramOptions::parse(int argc, char** argv) -> std::optional<ProgramOptions> {
    ProgramOptions result {};

    for (int opt; (opt = getopt(argc, argv, "m:i:b:w:")) != -1;) {
        switch (opt) {
        case 'm':
            result.matrix_size = std::atoi(optarg);
            break;

        case 'i':
            result.iterations = std::atoi(optarg);
            break;

        case 'b':
            result.block_size = std::atoi(optarg);
            break;

        case 'w':
            result.neighbor_width = std::atoi(optarg);
            break;

        case ':':
            std::fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            return std::nullopt;

        case '?':
            if (std::isprint(optopt))
                std::fprintf(stderr, "Unknown option `-%c'.\n", optopt);
            else
                std::fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            return std::nullopt;
        }
    }

    // Check the option arguments.
    if (result.matrix_size == 0) {
        std::fprintf(stderr, "Invalid matrix side length: %u\n", result.matrix_size);
        return std::nullopt;
    }

    if (result.block_size == 0) {
        std::fprintf(stderr, "Invalid block size: %u\n", result.block_size);
        return std::nullopt;
    }

    if (result.iterations == 0) {
        std::fprintf(stderr, "Invalid number of iterations: %u\n", result.iterations);
        return std::nullopt;
    }

    if (result.neighbor_width == 0 || result.neighbor_width * 2 + 1 > result.block_size) {
        std::fprintf(stderr, "Invalid width of neighboorhood: %u\n", result.neighbor_width);
    }

    return result;
}
