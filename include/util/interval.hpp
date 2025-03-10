#pragma once

#include <array>

namespace pmondriaan {

struct interval {
    long low;
    long high;
    long length() { return high - low; }
    long operator()(long part) {
        if (part == 0) {
            return low;
        } else {
            return high;
        }
    }
};

} // namespace pmondriaan
