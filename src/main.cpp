#include <athread.h>
#include <chrono>
#include <iostream>
#include <ratio>

#include "stencil/program_options.hpp"
#include "stencil/stencil.hpp"

int main(int argc, char** argv) {
    athread_init();

    // Parse program arugments.
    if (auto options = ProgramOptions::parse(argc, argv)) {
        constexpr unsigned count = 5;
        std::chrono::steady_clock::duration dma_duration = {};

        for (unsigned i = 0; i != count; ++i) {
            Stencil stencil(*options);
            auto const duration = stencil.run(Stencil::DMA);
            dma_duration += duration;
            std::cout << "DMA Method spent "
                      << static_cast<std::chrono::duration<double, std::milli>>(duration).count()
                      << "ms for " << options->iterations << " iterations.\n";
        }
        std::cout
            << "The average time taken by DMA method is "
            << static_cast<std::chrono::duration<double, std::milli>>(dma_duration / count).count()
            << "ms for " << options->iterations << " iterations.\n";

#if 0
        std::chrono::steady_clock::duration rma_duration = {};
        for (unsigned i = 0; i != count; ++i) {
            Stencil stencil(*options);
            auto const duration = stencil.run(Stencil::RMA);
            rma_duration += duration;
            std::cout << "RMA Method spent "
                      << static_cast<std::chrono::duration<double, std::milli>>(duration).count()
                      << "ms for " << options->iterations << " iterations.\n";
        }
        std::cout
            << "The average time taken by RMA method is "
            << static_cast<std::chrono::duration<double, std::milli>>(rma_duration / count).count()
            << "ms for " << options->iterations << " iterations.\n";
#endif
    } else {
        return 1;
    }
}
