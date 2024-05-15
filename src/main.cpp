#include <chrono>
#include <iostream>
#include <ratio>
#include <string_view>

#include <stencil/program_options.hpp>
#include <stencil/stencil.hpp>

#include <athread.h>

namespace {
void run_test(std::string_view method_name, ProgramOptions const& options) {
    std::chrono::steady_clock::duration total_duration = {};

    for (unsigned i = 0; i != options.repeat_count; ++i) {
        Stencil stencil(options);
        auto const duration = stencil.run(method_name);
        total_duration += duration;
        std::cout << method_name << " Method spent "
                  << static_cast<std::chrono::duration<double, std::milli>>(duration).count()
                  << "ms for " << options.iterations << " iterations.\n";
    }

    std::cout << "The average time taken by " << method_name << " method is "
              << static_cast<std::chrono::duration<double, std::milli>>(
                     total_duration / options.repeat_count
                 )
                     .count()
              << "ms for " << options.iterations << " iterations.\n";
}

void run_all_test(ProgramOptions const& options) {
    for (std::string const& method_name : options.method_names) {
        run_test(method_name, options);
    }
}
}  // namespace

int main(int argc, char** argv) {
    athread_init();

    // Parse program arugments.
    if (auto options = ProgramOptions::parse(argc, argv)) {
        run_all_test(*options);
    } else {
        return 1;
    }
}
