#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace VofFlow {
    template<typename T>
    std::size_t findClosestIndex(const T value, const std::vector<T>& list) {
        if (list.empty()) {
            throw std::runtime_error("List is empty!");
        }
        std::size_t idx = 0;
        T minDist = std::abs(list[0] - value);
        for (std::size_t i = 1; i < list.size(); i++) {
            const T dist = std::abs(list[i] - value);
            if (dist < minDist) {
                minDist = dist;
                idx = i;
            }
        }
        return idx;
    }
} // namespace VofFlow
