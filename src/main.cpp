#include <athread.h>
#include <chrono>
#include <iostream>
#include <ratio>
#include <string_view>

#include "stencil/program_options.hpp"
#include "stencil/stencil.hpp"

namespace {
void run_test(
    Stencil::InputMethod method,
    std::string_view name,
    unsigned count,
    ProgramOptions const& options
) {
    std::chrono::steady_clock::duration total_duration = {};

    for (unsigned i = 0; i != count; ++i) {
        Stencil stencil(options);
        auto const duration = stencil.run(method);
        total_duration += duration;
        std::cout << name << " Method spent "
                  << static_cast<std::chrono::duration<double, std::milli>>(duration).count()
                  << "ms for " << options.iterations << " iterations.\n";
    }

    std::cout
        << "The average time taken by " << name << " method is "
        << static_cast<std::chrono::duration<double, std::milli>>(total_duration / count).count()
        << "ms for " << options.iterations << " iterations.\n";
}
}  // namespace

int main(int argc, char** argv) {
    athread_init();

    // Parse program arugments.
    if (auto options = ProgramOptions::parse(argc, argv)) {
        constexpr unsigned count = 5;

        run_test(Stencil::DMA, "DMA", count, *options);
        run_test(Stencil::DMA_STATIC_UNROLL, "DMA Static Unroll", count, *options);
        run_test(Stencil::DMA_SLAVE_PACK, "DMA Slave Pack", count, *options);
        run_test(Stencil::RMA, "RMA", count, *options);
    } else {
        return 1;
    }
}
