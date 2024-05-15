#include <optional>

#include <stencil/program_options.hpp>

#include <CLI/CLI.hpp>
#include <CLI/Error.hpp>

auto ProgramOptions::parse(int argc, char** argv) -> std::optional<ProgramOptions> {
    CLI::App app;

    ProgramOptions result {};

    app.add_option("-s,--matrix-size", result.matrix_size, "The side length of the input matrix.")
        ->required();
    app.add_option("-i,--iteration", result.iterations, "The number of iterations.")->required();
    app.add_option(
           "-b,--block-size",
           result.block_size,
           "The side length of the block into which the matrix is divided."
    )
        ->required();
    app.add_option("-r,--radius", result.radius, "The radius of the stencil shape.")
        ->default_val(1);
    app.add_option(
           "-R,--repeat",
           result.repeat_count,
           "The number of replicate runs for each method."
    )
        ->default_val(1);
    app.add_option("-m,--methods", result.method_names, "List of methods to be tested.")
        ->required();

    try {
        app.parse(argc, argv);
    } catch (CLI::ParseError const& e) {
        app.exit(e);
        return std::nullopt;
    }

    return result;
}
