#ifndef STENCIL_PROGRAM_OPTIONS_HPP
#define STENCIL_PROGRAM_OPTIONS_HPP

#include <optional>
#include <string>
#include <vector>

struct ProgramOptions {
    /// Side length of the matrix.
    unsigned matrix_size;
    /// The side length of the block obtained by each CPE. If 64 CPEs are not enough to process the
    /// entire matrix in one round, then they divide the matrix into large chunks and process them
    /// in multiple rounds.
    unsigned block_size;
    /// Number of iterations.
    unsigned iterations;
    /// The radius of the stencil shape.
    unsigned radius;
    /// Number of replicate runs for each method.
    unsigned repeat_count;
    /// The list of names of methods to be tested.
    std::vector<std::string> method_names;

    /// Parse the program arguments.
    static auto parse(int argc, char** argv) -> std::optional<ProgramOptions>;
};

#endif  // STENCIL_PROGRAM_OPTIONS_HPP
