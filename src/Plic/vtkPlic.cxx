#include "vtkPlic.h"

#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>

#include "Grid/DomainInfo.h"
#include "Grid/GridIterator.h"
#include "Grid/GridTypes.h"
#include "Misc/Profiling.h"
#include "Plic/Plic.h"
#include "Plic/PlicPolyData.h"

vtkStandardNewMacro(vtkPlic);

int vtkPlic::calcPlic(vtkRectilinearGrid* input, vtkPolyData* output) {
    ZoneScoped;

    VofFlow::DomainInfo domainInfo(input, mpiController_);

    vtkDataArray* vofData = GetInputArrayToProcess(0, input);

    VofFlow::vtkPlicPolyData vtkPlic;

    for (const VofFlow::gridCoords_t& e_coords : VofFlow::GridRange(domainInfo.localZeroExtent())) {
        const double f = vofData->GetComponent(domainInfo.localExtentCoordToIdx(e_coords), 0);
        if (f < Epsilon || f > 1.0 - Epsilon) {
            continue;
        }

        const auto& plic = VofFlow::calcPlicPoly(domainInfo, domainInfo.localExtentCoordToGridCoord(e_coords), vofData,
            Epsilon, NumIterations);
        vtkPlic.addPolygon(plic);
    }

    output->ShallowCopy(vtkPlic.getPolyData());

    return 1;
}
