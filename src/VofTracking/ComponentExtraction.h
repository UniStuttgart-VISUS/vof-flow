#pragma once

#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

#include <vtkDataArray.h>
#include <vtkIntArray.h>
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>

#include "Grid/DomainInfo.h"
#include "Grid/GridTypes.h"

namespace VofFlow {
    struct ComponentResult {
        vtkSmartPointer<vtkIntArray> gridLabels = nullptr;
        int numLabels = 0;
    };

    class ComponentExtractor {
    public:
        explicit ComponentExtractor(const DomainInfo& domainInfo, vtkMPIController* mpiController = nullptr)
            : domainInfo_(domainInfo),
              mpiController_(mpiController) {}

        ComponentResult extractComponents(const vtkSmartPointer<vtkDataArray>& vof);

        ComponentResult checkComponents(const vtkSmartPointer<vtkDataArray>& inComponents);

    private:
        static inline uint64_t makeCombinedLabel(int processId, int label) {
            return (static_cast<uint64_t>(processId) << 32) + static_cast<uint32_t>(label);
        }

        static inline int getLabelProcessId(uint64_t l) {
            return static_cast<int>(static_cast<uint32_t>(l >> 32));
        }

        static inline int getLabelIdx(uint64_t l) {
            return static_cast<int>(static_cast<uint32_t>(l & 0xFFFFFFFF));
        }

        static inline std::tuple<int, int> splitCombinedLabel(uint64_t l) {
            return {
                getLabelProcessId(l),
                getLabelIdx(l),
            };
        }

        struct BorderExtent {
            int neighborId;
            extent_t extent;

            inline bool operator<(const BorderExtent& other) const {
                return neighborId < other.neighborId;
            }
        };

        struct LabelPair {
            uint64_t label1;
            uint64_t label2;

            LabelPair() = default;
            LabelPair(uint64_t l1, uint64_t l2) : label1(l1), label2(l2) {
                // Label pairs should have smaller label first to avoid duplicates.
                if (label2 < label1) {
                    std::swap(label1, label2);
                }
            }

            inline bool operator<(const LabelPair& other) const {
                return (label1 == other.label1) ? (label2 < other.label2) : (label1 < other.label1);
            }
        };

        struct MappingPair {
            int from;
            int to;
        };

        ComponentResult extractComponentsLocal(const vtkSmartPointer<vtkDataArray>& vof);

        std::tuple<std::vector<BorderExtent>, std::vector<BorderExtent>> calcBorderExtents();

        const DomainInfo& domainInfo_;
        vtkMPIController* mpiController_;
    };
} // namespace VofFlow
