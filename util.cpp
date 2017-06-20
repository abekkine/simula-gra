#include "util.h"

#include <iostream>

namespace util {

    auto t0 = std::chrono::steady_clock::now();

    void start_clock() {
        t0 = std::chrono::steady_clock::now();
    }

    void stop_clock() {
        auto t1 = std::chrono::steady_clock::now();
        auto diff = t1 - t0;
#ifdef IN_SECONDS
        std::cout << std::chrono::duration <double> (diff).count();
        std::cout << " s" << std::endl;
#else
        std::cout << std::chrono::duration <double, std::milli> (diff).count();
        std::cout << " ms" << std::endl;
#endif
    }
}

