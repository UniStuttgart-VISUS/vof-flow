#pragma once

#include <algorithm>

#include "GridTypes.h"

namespace VofFlow {
    class GridRange {
    public:
        GridRange(const gridCoords_t& startPos, const gridCoords_t& endPos) : startPos_(startPos), endPos_(endPos) {}

        explicit GridRange(const extent_t& extent) {
            startPos_ = {extent[0], extent[2], extent[4]};
            endPos_ = {extent[1], extent[3], extent[5]};
        }

        GridRange(const dim_t& dimensions, const gridCoords_t& start, int searchRange) {
            startPos_ = {
                std::max(start[0] - searchRange, 0),
                std::max(start[1] - searchRange, 0),
                std::max(start[2] - searchRange, 0),
            };
            endPos_ = {
                std::min(start[0] + searchRange + 1, dimensions[0]),
                std::min(start[1] + searchRange + 1, dimensions[1]),
                std::min(start[2] + searchRange + 1, dimensions[2]),
            };
        }

        class Iterator {
        public:
            Iterator(const gridCoords_t& startPos, const gridCoords_t& endPos)
                : position_(startPos),
                  startPos_(startPos),
                  endPos_(endPos) {
                if (position_[0] >= endPos_[0] || position_[1] >= endPos_[1] || position_[2] >= endPos_[2]) {
                    position_ = endPos_;
                }
            }

            inline Iterator& operator++() {
                position_[0]++;
                if (position_[0] >= endPos_[0]) {
                    position_[0] = startPos_[0];
                    position_[1]++;
                    if (position_[1] >= endPos_[1]) {
                        position_[1] = startPos_[1];
                        position_[2]++;
                        if (position_[2] >= endPos_[2]) {
                            position_ = endPos_;
                        }
                    }
                }
                return *this;
            }

            inline const gridCoords_t& operator*() const {
                return position_;
            }

            inline bool operator==(const Iterator& other) {
                return position_ == other.position_;
            }

            inline bool operator!=(const Iterator& other) {
                return position_ != other.position_;
            }

        private:
            gridCoords_t position_;
            const gridCoords_t& startPos_;
            const gridCoords_t& endPos_;
        };

        [[nodiscard]] inline auto begin() const {
            return Iterator(startPos_, endPos_);
        }

        [[nodiscard]] inline auto end() const {
            return Iterator(endPos_, endPos_);
        }

    private:
        gridCoords_t startPos_;
        gridCoords_t endPos_;
    };

} // namespace VofFlow
