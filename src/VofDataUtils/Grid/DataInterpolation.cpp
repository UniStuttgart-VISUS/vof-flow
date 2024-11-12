#include "DataInterpolation.h"

#include <algorithm>
#include <array>
#include <stdexcept>

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkAssume.h>
#include <vtkDataArrayAccessor.h>
#include <vtkType.h>

namespace {
    struct CellDataInterpolation3Worker {

        const VofFlow::DomainInfo& domainInfo_;
        const VofFlow::gridCoords_t& cellCoords_;
        const VofFlow::posCoords_t& pCoords_;
        VofFlow::vec3 result_;

        CellDataInterpolation3Worker(const VofFlow::DomainInfo& domainInfo, const VofFlow::gridCoords_t& cellCoords,
            const VofFlow::posCoords_t& pCoords)
            : domainInfo_(domainInfo),
              cellCoords_(cellCoords),
              pCoords_(pCoords),
              result_(0.0f, 0.0f, 0.0f) {}

        template<typename ValueType>
        void operator()(vtkAOSDataArrayTemplate<ValueType>* dataArray) {
            VTK_ASSUME(dataArray->GetNumberOfComponents() == 3);

            vtkDataArrayAccessor<vtkAOSDataArrayTemplate<ValueType>> data(dataArray);

            int lx = cellCoords_[0];
            int ly = cellCoords_[1];
            int lz = cellCoords_[2];

            double x = pCoords_[0] - 0.5;
            double y = pCoords_[1] - 0.5;
            double z = pCoords_[2] - 0.5;

            if (pCoords_[0] < 0.5) {
                lx -= 1;
                x = pCoords_[0] + 0.5;
            }
            if (pCoords_[1] < 0.5) {
                ly -= 1;
                y = pCoords_[1] + 0.5;
            }
            if (pCoords_[2] < 0.5) {
                lz -= 1;
                z = pCoords_[2] + 0.5;
            }

            int ux = lx + 1;
            int uy = ly + 1;
            int uz = lz + 1;

            lx = std::max(lx, 0);
            ly = std::max(ly, 0);
            lz = std::max(lz, 0);
            ux = std::min(ux, domainInfo_.cellDims()[0] - 1);
            uy = std::min(uy, domainInfo_.cellDims()[1] - 1);
            uz = std::min(uz, domainInfo_.cellDims()[2] - 1);

            const std::array<vtkIdType, 8> idx{
                domainInfo_.gridCoordToIdx(lx, ly, lz),
                domainInfo_.gridCoordToIdx(ux, ly, lz),
                domainInfo_.gridCoordToIdx(lx, uy, lz),
                domainInfo_.gridCoordToIdx(ux, uy, lz),
                domainInfo_.gridCoordToIdx(lx, ly, uz),
                domainInfo_.gridCoordToIdx(ux, ly, uz),
                domainInfo_.gridCoordToIdx(lx, uy, uz),
                domainInfo_.gridCoordToIdx(ux, uy, uz),
            };
            std::array<VofFlow::vec3, 8> values;
            for (int i = 0; i < 8; i++) {
                values[i] = VofFlow::vec3{
                    static_cast<float>(data.Get(idx[i], 0)),
                    static_cast<float>(data.Get(idx[i], 1)),
                    static_cast<float>(data.Get(idx[i], 2)),
                };
            }

            const VofFlow::vec3 a = VofFlow::lerp(values[0], values[1], static_cast<float>(x));
            const VofFlow::vec3 b = VofFlow::lerp(values[2], values[3], static_cast<float>(x));
            const VofFlow::vec3 c = VofFlow::lerp(values[4], values[5], static_cast<float>(x));
            const VofFlow::vec3 d = VofFlow::lerp(values[6], values[7], static_cast<float>(x));

            const VofFlow::vec3 e = VofFlow::lerp(a, b, static_cast<float>(y));
            const VofFlow::vec3 f = VofFlow::lerp(c, d, static_cast<float>(y));

            result_ = VofFlow::lerp(e, f, static_cast<float>(z));
        }
    };
} // namespace

VofFlow::vec3 VofFlow::interpolateCellDataVec3(const vtkSmartPointer<vtkDataArray>& data, const DomainInfo& domainInfo,
    const gridCoords_t& cellCoords, const posCoords_t& pCoords) {
    CellDataInterpolation3Worker worker(domainInfo, cellCoords, pCoords);

    typedef vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals> Dispatcher;
    if (!Dispatcher::Execute(data, worker)) {
        throw std::runtime_error("Cannot dispatch array worker!");
    }

    return worker.result_;
}
