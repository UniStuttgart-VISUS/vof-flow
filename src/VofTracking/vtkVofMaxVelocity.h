#pragma once

#include <memory>

#include <vtkMPIController.h>
#include <vtkPVVersion.h>
#include <vtkTableAlgorithm.h>
#include <vtkTemporalAlgorithm.h>
#include <vtkType.h>

#include "Math/Vector.h"

#include "VofTrackingModule.h"

#ifndef __VTK_WRAP__
#define vtkTableAlgorithm vtkTemporalAlgorithm<vtkTableAlgorithm>
#endif

struct MaxResult;
struct HistResult;

class VOFTRACKING_EXPORT vtkVofMaxVelocity : public vtkTableAlgorithm {
public:
    vtkTypeMacro(vtkVofMaxVelocity, vtkTableAlgorithm);
#ifndef __VTK_WRAP__
#undef vtkTableAlgorithm
#endif
    static vtkVofMaxVelocity* New();
    void PrintSelf(ostream& os, vtkIndent indent) override;

#if defined(__VTK_WRAP__) || defined(__WRAP_GCCXML)
    vtkCreateWrappedTemporalAlgorithmInterface();
#endif

    // ParaView GUI
    vtkGetMacro(Epsilon, double);
    vtkSetMacro(Epsilon, double);

    vtkVofMaxVelocity(const vtkVofMaxVelocity&) = delete;
    void operator=(const vtkVofMaxVelocity&) = delete;

protected:
    vtkVofMaxVelocity();
    ~vtkVofMaxVelocity() override = default;

    int FillInputPortInformation(int port, vtkInformation* info) override;
#if PARAVIEW_VERSION_NUMBER < PARAVIEW_VERSION_CHECK(5, 13, 0)
    int RequestInformation(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;
#endif

    int Initialize(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;
    int Execute(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;
    int Finalize(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;

private:
    [[nodiscard]] inline bool isRank0() const {
        return mpiController_ == nullptr || mpiController_->GetLocalProcessId() == 0;
    }

    // ParaView GUI
    double Epsilon;

    // Internal State
    std::unique_ptr<MaxResult> maxResult_;
    std::unique_ptr<HistResult> histResult_;
    VofFlow::vec3 cellDimsUniform_;
    vtkMPIController* mpiController_;

    static constexpr vtkIdType NUM_BINS = 128;

#if PARAVIEW_VERSION_NUMBER < PARAVIEW_VERSION_CHECK(5, 13, 0)
    // Parent InputTimeSteps is a private member in ParaView 5.12. We need to get our own copy from input data.
    std::vector<double> InputTimeSteps;
#endif
};
