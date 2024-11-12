#include "vtkPlic3.h"

#include <vector>

#include <CGAL/intersections.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>

#include "Grid/DomainInfo.h"
#include "Grid/GridIterator.h"
#include "Grid/GridTypes.h"
#include "Misc/CgalUtil.h"
#include "Misc/Profiling.h"
#include "Misc/VofData.h"
#include "Plic/Plic3.h"
#include "Plic/PlicPolyData.h"
#include "Plic/Polygon3D.h"

vtkStandardNewMacro(vtkPlic3);

int vtkPlic3::calcPlic(vtkRectilinearGrid* input, vtkPolyData* output) {
    ZoneScoped;

    VofFlow::DomainInfo domainInfo(input, mpiController_);

    VofFlow::VofData vofData{
        GetInputArrayToProcess(0, input),
        GetInputArrayToProcess(1, input),
        GetInputArrayToProcess(2, input),
    };

    VofFlow::vtkPlic3PolyData vtkPlic;

    for (const VofFlow::gridCoords_t& e_coords : VofFlow::GridRange(domainInfo.localZeroExtent())) {
        const auto& plic = VofFlow::calcPlic3CellClass(domainInfo, domainInfo.localExtentCoordToGridCoord(e_coords),
            vofData, Epsilon, NumIterations);

        if (plic.cellClass == VofFlow::CellClass::EMPTY || plic.cellClass == VofFlow::CellClass::FULL_PHASE1 ||
            plic.cellClass == VofFlow::CellClass::FULL_PHASE2) {
            continue;
        }

        // PhaseType: If interface of first phase value is 1, otherwise if interface of second phase, value is 2.
        // CellClass must be one of: INTERFACE_PH1_PH2, INTERFACE_PHASE1, INTERFACE_PHASE2, INTERFACE_ALL
        const unsigned char phaseType = (plic.cellClass == VofFlow::CellClass::INTERFACE_PHASE2) ? 2 : 1;

        const auto& plic1 = plic.plic1.value();
        vtkPlic.addPolygon(plic.cellClass, phaseType, {plic1.cell.getPolygon(0), plic1.iter, plic1.err});

        if (plic.cellClass == VofFlow::CellClass::INTERFACE_ALL) {
            const auto& plic2 = plic.plic2.value();
            const auto& plane1 = plic1.cell.getPlane(0);
            const auto& plane2 = plic2.cell.getPlane(1);

            std::vector<VofFlow::K::Point_3> polyPoints;
            const auto& poly = plic2.cell.getPolygon(1); // Store temporary object
            for (const auto& p : poly.points()) {
                if (plane1.has_on_positive_side(p)) {
                    polyPoints.push_back(p);
                }
            }

            const auto result = CGAL::intersection(plane1, plane2);
            if (result) {
                if (const VofFlow::K::Line_3* l = boost::get<VofFlow::K::Line_3>(&*result)) {
                    const auto result2 = CGAL::intersection(plic1.cell.cellCube(), *l);
                    if (result2) {
                        if (const VofFlow::K::Segment_3* s = boost::get<VofFlow::K::Segment_3>(&*result2)) {
                            polyPoints.push_back(s->source());
                            polyPoints.push_back(s->target());
                        } else {
                            const VofFlow::K::Point_3* p = boost::get<VofFlow::K::Point_3>(&*result2);
                            polyPoints.push_back(*p);
                        }
                    }
                } else {
                    // Planes are identical. No need to add a second plane at the very same position.
                }
            } else {
                // Parallel planes. `has_on_positive_side` has either added all or none points.
            }

            if (polyPoints.size() >= 3) {
                // Second plane in three-phase cell has always PhaseType 2
                vtkPlic.addPolygon(plic.cellClass, 2, {VofFlow::Polygon3D(plane2, polyPoints), plic2.iter, plic2.err});
            }
        }
    }

    output->ShallowCopy(vtkPlic.getPolyData());

    return 1;
}
