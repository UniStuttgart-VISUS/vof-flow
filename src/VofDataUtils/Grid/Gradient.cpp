#include "Gradient.h"

#include <stdexcept>

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkAssume.h>
#include <vtkDataArrayAccessor.h>

namespace {
    struct GradientWorker {
        const VofFlow::DomainInfo& domainInfo_;
        const int x_;
        const int y_;
        const int z_;
        VofFlow::vec3 gradient_;
        const bool borderXmin_;
        const bool borderXmax_;
        const bool borderYmin_;
        const bool borderYmax_;
        const bool borderZmin_;
        const bool borderZmax_;

        GradientWorker(const VofFlow::DomainInfo& domainInfo, const VofFlow::gridCoords_t& g_coords)
            : domainInfo_(domainInfo),
              x_(g_coords[0]),
              y_(g_coords[1]),
              z_(g_coords[2]),
              gradient_(0.0f, 0.0f, 0.0f),
              borderXmin_(x_ <= 0),
              borderXmax_(x_ >= domainInfo_.cellDims()[0] - 1),
              borderYmin_(y_ <= 0),
              borderYmax_(y_ >= domainInfo_.cellDims()[1] - 1),
              borderZmin_(z_ <= 0),
              borderZmax_(z_ >= domainInfo_.cellDims()[2] - 1) {}

        template<typename ValueType>
        inline ValueType gradientXDim(const vtkDataArrayAccessor<vtkAOSDataArrayTemplate<ValueType>>& data, int y,
            int z) {
            const auto& cellSizesX = domainInfo_.cellSizesX();
            if (borderXmin_) {
                // forward difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x_, y, z), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x_ + 1, y, z), 0);
                return (val1 - val0) / (0.5 * (cellSizesX[x_] + cellSizesX[x_ + 1]));
            } else if (borderXmax_) {
                // backward difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x_ - 1, y, z), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x_, y, z), 0);
                return (val1 - val0) / (0.5 * (cellSizesX[x_ - 1] + cellSizesX[x_]));
            } else {
                // central difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x_ - 1, y, z), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x_ + 1, y, z), 0);
                return (val1 - val0) / (0.5 * cellSizesX[x_ - 1] + cellSizesX[x_] + 0.5 * cellSizesX[x_ + 1]);
            }
        }

        template<typename ValueType>
        inline ValueType gradientYDim(const vtkDataArrayAccessor<vtkAOSDataArrayTemplate<ValueType>>& data, int x,
            int z) {
            const auto& cellSizesY = domainInfo_.cellSizesY();
            if (borderYmin_) {
                // forward difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x, y_, z), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x, y_ + 1, z), 0);
                return (val1 - val0) / (0.5 * (cellSizesY[y_] + cellSizesY[y_ + 1]));
            } else if (borderYmax_) {
                // backward difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x, y_ - 1, z), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x, y_, z), 0);
                return (val1 - val0) / (0.5 * (cellSizesY[y_ - 1] + cellSizesY[y_]));
            } else {
                // central difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x, y_ - 1, z), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x, y_ + 1, z), 0);
                return (val1 - val0) / (0.5 * cellSizesY[y_ - 1] + cellSizesY[y_] + 0.5 * cellSizesY[y_ + 1]);
            }
        }

        template<typename ValueType>
        inline ValueType gradientZDim(const vtkDataArrayAccessor<vtkAOSDataArrayTemplate<ValueType>>& data, int x,
            int y) {
            const auto& cellSizesZ = domainInfo_.cellSizesZ();
            if (borderZmin_) {
                // forward difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x, y, z_), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x, y, z_ + 1), 0);
                return (val1 - val0) / (0.5 * (cellSizesZ[z_] + cellSizesZ[z_ + 1]));
            } else if (borderZmax_) {
                // backward difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x, y, z_ - 1), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x, y, z_), 0);
                return (val1 - val0) / (0.5 * (cellSizesZ[z_ - 1] + cellSizesZ[z_]));
            } else {
                // central difference
                const ValueType val0 = data.Get(domainInfo_.gridCoordToIdx(x, y, z_ - 1), 0);
                const ValueType val1 = data.Get(domainInfo_.gridCoordToIdx(x, y, z_ + 1), 0);
                return (val1 - val0) / (0.5 * cellSizesZ[z_ - 1] + cellSizesZ[z_] + 0.5 * cellSizesZ[z_ + 1]);
            }
        }

        template<typename ValueType>
        void operator()(vtkAOSDataArrayTemplate<ValueType>* dataArray) {
            VTK_ASSUME(dataArray->GetNumberOfComponents() == 1);

            vtkDataArrayAccessor<vtkAOSDataArrayTemplate<ValueType>> data(dataArray);

            // TODO weighting by cell size is not implemented yet!

            // X
            ValueType gradientX = 0.0;
            {
                ValueType weightsX = 0.0;
                for (int y = 0; y < 3; y++) {
                    if ((y == 0 && borderYmin_) || (y == 2 && borderYmax_)) {
                        continue;
                    }
                    const int weightY = (y == 1) ? 2 : 1;
                    for (int z = 0; z < 3; z++) {
                        if ((z == 0 && borderZmin_) || (z == 2 && borderZmax_)) {
                            continue;
                        }
                        const int weightZ = (z == 1) ? 2 : 1;
                        const double weight = static_cast<double>(weightY * weightZ);

                        gradientX += weight * gradientXDim(data, y_ - 1 + y, z_ - 1 + z);
                        weightsX += weight;
                    }
                }
                gradientX /= weightsX;
            }

            // Y
            ValueType gradientY = 0.0;
            {
                ValueType weightsY = 0.0;
                for (int x = 0; x < 3; x++) {
                    if ((x == 0 && borderXmin_) || (x == 2 && borderXmax_)) {
                        continue;
                    }
                    const int weightX = (x == 1) ? 2 : 1;
                    for (int z = 0; z < 3; z++) {
                        if ((z == 0 && borderZmin_) || (z == 2 && borderZmax_)) {
                            continue;
                        }
                        const int weightZ = (z == 1) ? 2 : 1;
                        const double weight = static_cast<double>(weightX * weightZ);

                        gradientY += weight * gradientYDim(data, x_ - 1 + x, z_ - 1 + z);
                        weightsY += weight;
                    }
                }
                gradientY /= weightsY;
            }

            // Z
            ValueType gradientZ = 0.0;
            {
                ValueType weightsZ = 0.0;
                for (int x = 0; x < 3; x++) {
                    if ((x == 0 && borderXmin_) || (x == 2 && borderXmax_)) {
                        continue;
                    }
                    const int weightX = (x == 1) ? 2 : 1;
                    for (int y = 0; y < 3; y++) {
                        if ((y == 0 && borderYmin_) || (y == 2 && borderYmax_)) {
                            continue;
                        }
                        const int weightY = (y == 1) ? 2 : 1;
                        const double weight = static_cast<double>(weightX * weightY);

                        gradientZ += weight * gradientZDim(data, x_ - 1 + x, y_ - 1 + y);
                        weightsZ += weight;
                    }
                }
                gradientZ /= weightsZ;
            }

            gradient_ = VofFlow::vec3{
                static_cast<float>(gradientX),
                static_cast<float>(gradientY),
                static_cast<float>(gradientZ),
            };
        }
    };
} // namespace

VofFlow::vec3 VofFlow::gradient(const DomainInfo& domainInfo, const gridCoords_t& g_coords,
    const vtkSmartPointer<vtkDataArray>& data) {
    GradientWorker worker(domainInfo, g_coords);

    typedef vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals> Dispatcher;
    if (!Dispatcher::Execute(data, worker)) {
        throw std::runtime_error("Cannot dispatch array worker!");
    }

    return worker.gradient_;
}
