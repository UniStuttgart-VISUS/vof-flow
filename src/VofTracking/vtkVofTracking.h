#pragma once

#include <memory>
#include <string>
#include <vector>

#include <vtkDataSetAlgorithm.h>
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>
#include <vtkTimeStamp.h>
#include <vtkType.h>

#include "VofTrackingModule.h"

class vtkDataArray;
class vtkImplicitFunction;

namespace VofFlow {
    class DomainInfo;
    class SeedCoordInfo;
    struct SeedPoints;
    struct Particles;
    struct OutOfBoundsParticles;
} // namespace VofFlow

#define vtkSetMacroNoInternalModified(name, type) \
    virtual void Set##name(type _arg) {           \
        if (this->name != _arg) {                 \
            this->name = _arg;                    \
            this->ModifiedNotInternally();        \
        }                                         \
    }

class VOFTRACKING_EXPORT vtkVofTracking : public vtkDataSetAlgorithm {
public:
    static vtkVofTracking* New();
    vtkTypeMacro(vtkVofTracking, vtkDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent) override;

    // ParaView GUI
    vtkGetMacro(UseThreePhase, int);
    vtkSetMacro(UseThreePhase, int);

    vtkGetMacro(UseComponents, int);
    vtkSetMacro(UseComponents, int);

    vtkGetMacro(UseTargetTimeStep, int);
    vtkSetMacro(UseTargetTimeStep, int);

    vtkGetMacro(InitTimeStep, int);
    vtkSetMacro(InitTimeStep, int);

    vtkGetMacro(TargetTimeStep, int);
    vtkSetMacro(TargetTimeStep, int);

    vtkGetMacro(Refinement, int);
    vtkSetMacro(Refinement, int);

    vtkGetMacro(NeighborCorrection, int);
    vtkSetMacro(NeighborCorrection, int);

    vtkGetMacro(CellCorrection, int);
    vtkSetMacro(CellCorrection, int);

    vtkGetMacro(PLICCorrection, int);
    vtkSetMacro(PLICCorrection, int);

    vtkGetMacro(IntegrationMethod, int);
    vtkSetMacro(IntegrationMethod, int);

    vtkGetMacro(IntegrationSubSteps, int);
    vtkSetMacro(IntegrationSubSteps, int);

    vtkGetMacro(Epsilon, double);
    vtkSetMacro(Epsilon, double);

    vtkGetMacro(NumIterations, int);
    vtkSetMacro(NumIterations, int);

    vtkGetMacro(GhostCells, int);
    vtkSetMacro(GhostCells, int);

    vtkGetMacro(CutLabels, int);
    vtkSetMacroNoInternalModified(CutLabels, int);

    vtkGetMacro(CutFunction, vtkImplicitFunction*);
    vtkSetMacroNoInternalModified(CutFunction, vtkImplicitFunction*);

    vtkGetMacro(BoundaryMethod, int);
    vtkSetMacro(BoundaryMethod, int);

    vtkGetMacro(OutputDataType, int);
    vtkSetMacro(OutputDataType, int);

    vtkGetMacro(OutputState, int);
    vtkSetMacro(OutputState, int);

    vtkGetMacro(OutputTimeMeasure, int);
    vtkSetMacro(OutputTimeMeasure, int);

    vtkGetMacro(MirrorXMin, int);
    vtkSetMacro(MirrorXMin, int);

    vtkGetMacro(MirrorXMax, int);
    vtkSetMacro(MirrorXMax, int);

    vtkGetMacro(MirrorYMin, int);
    vtkSetMacro(MirrorYMin, int);

    vtkGetMacro(MirrorYMax, int);
    vtkSetMacro(MirrorYMax, int);

    vtkGetMacro(MirrorZMin, int);
    vtkSetMacro(MirrorZMin, int);

    vtkGetMacro(MirrorZMax, int);
    vtkSetMacro(MirrorZMax, int);

    // Override GetMTime to also check for CutFunction changes.
    vtkMTimeType GetMTime() override;

    // We want to be able to change some of our filter properties without needing to rerun our full algorithm.
    // Therefore, we use internalMTime_ to track state in addition the MTime of this vtkObject. To be notified of all
    // external changes we need to override Modified. In addition, ModifiedNotInternally only updates the MTime of the
    // vtkObject to let ParaView know there was a change, but without updating the internal MTime to not rerun the
    // algorithm.
    void Modified() override;
    void ModifiedNotInternally();

    vtkVofTracking(const vtkVofTracking&) = delete;
    void operator=(const vtkVofTracking&) = delete;

protected:
    vtkVofTracking();
    ~vtkVofTracking() override;

    int FillInputPortInformation(int port, vtkInformation* info) override;
    int FillOutputPortInformation(int port, vtkInformation* info) override;
    int RequestInformation(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;
    int RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;
    int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;
    int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
        vtkInformationVector* outputVector) override;

private:
    [[nodiscard]] inline bool hasMPI() const {
        return mpiController_ != nullptr && mpiController_->GetCommunicator() != nullptr;
    }

    [[nodiscard]] inline bool isRank0() const {
        return mpiController_ == nullptr || mpiController_->GetLocalProcessId() == 0;
    }

    [[nodiscard]] std::string stateJson() const;

    enum State {
        NONE = 0,
        INIT,
        ADVECTION,
        OUTPUT,
    };

    static inline std::string toString(State s) {
        switch (s) {
            case NONE:
                return "NONE";
            case INIT:
                return "INIT";
            case ADVECTION:
                return "ADVECTION";
            case OUTPUT:
                return "OUTPUT";
        }
        return "INVALID STATE";
    }

    // ParaView GUI
    int UseThreePhase;
    int UseComponents;
    int UseTargetTimeStep;
    int InitTimeStep;
    int TargetTimeStep;
    int Refinement;
    int NeighborCorrection;
    int CellCorrection;
    int PLICCorrection;
    int IntegrationMethod;
    int IntegrationSubSteps;
    double Epsilon;
    int NumIterations;
    int GhostCells;
    int CutLabels;
    vtkImplicitFunction* CutFunction;
    int BoundaryMethod;
    int OutputDataType;
    int OutputState;
    int OutputTimeMeasure;
    int MirrorXMin;
    int MirrorXMax;
    int MirrorYMin;
    int MirrorYMax;
    int MirrorZMin;
    int MirrorZMax;

    // Time steps + state
    std::vector<double> timeStepValues_;
    int timestepT0_;
    int timestepT1_;
    int timestepInit_;
    int timestepTarget_;
    State internalState_;
    vtkTimeStamp internalMTime_;
    vtkMTimeType lastMTime_;

    // MPI / domain context
    vtkMPIController* mpiController_;
    const int requiredNumGhostLevels_;
    std::unique_ptr<VofFlow::DomainInfo> domainInfo_;
    std::unique_ptr<VofFlow::SeedCoordInfo> seedInfo_;

    // Data from last time step
    vtkSmartPointer<vtkDataArray> lastDataVelocity_;

    // Particles
    std::unique_ptr<VofFlow::SeedPoints> seedPoints_;
    std::unique_ptr<VofFlow::Particles> particles_;
    std::unique_ptr<VofFlow::OutOfBoundsParticles> oobParticles_;
};
