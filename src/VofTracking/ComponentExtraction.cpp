#include "ComponentExtraction.h"

#include <algorithm>
#include <cstddef>
#include <queue>
#include <set>
#include <stack>
#include <stdexcept>
#include <unordered_map>

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>
#include <vtkAssume.h>
#include <vtkCommunicator.h>
#include <vtkDataArray.h>
#include <vtkDataArrayAccessor.h>
#include <vtkMPICommunicator.h>
#include <vtkType.h>

#include "Constants.h"
#include "Grid/GridData.h"
#include "Misc/Profiling.h"
#include "VtkUtil/VtkUtilArray.h"
#include "VtkUtil/VtkUtilMPI.h"

namespace {
    struct ExtractComponentsWorker {
        int nextLabelId_;
        double epsilon_;
        const VofFlow::DomainInfo& domainInfo_;
        vtkSmartPointer<vtkIntArray> labels_;
        std::vector<vtkIdType> labelMinCellIdx_;

        explicit ExtractComponentsWorker(const VofFlow::DomainInfo& domainInfo)
            : nextLabelId_(0),
              epsilon_(0.0),
              domainInfo_(domainInfo) {}

        [[nodiscard]] inline bool isGhost(vtkIdType idx) const {
            return domainInfo_.ghostArray() != nullptr && domainInfo_.ghostArray()->GetValue(idx) > 0;
        }

        template<typename ValueType>
        void operator()(vtkAOSDataArrayTemplate<ValueType>* vof) {
            ZoneScoped;

            if (vof->GetNumberOfComponents() != 1) {
                throw std::runtime_error("Invalid number of components in vof array!");
            }
            VTK_ASSUME(vof->GetNumberOfComponents() == 1);

            labels_ = VofFlow::createVtkArray<vtkIntArray>(VofFlow::ArrayNames::LABELS, 1, vof->GetNumberOfTuples(),
                VofFlow::ErrorLabels::EMPTY_COMPONENT);

            vtkDataArrayAccessor<vtkAOSDataArrayTemplate<ValueType>> vofData(vof);
            std::stack<vtkIdType> idxStack;
            const vtkIdType numCells = domainInfo_.numCells();
            vtkIdType idx = 0;

            const int offsetX = 1;
            const int offsetY = domainInfo_.cellDims()[0];
            const int offsetZ = domainInfo_.cellDims()[0] * domainInfo_.cellDims()[1];

            while (idx < numCells) {
                while (idx < numCells &&
                       (vofData.Get(idx, 0) <= epsilon_ || labels_->GetValue(idx) >= 0 || isGhost(idx))) {
                    idx++;
                }
                if (idx >= numCells) {
                    break;
                }
                // We loop over increasing idx, so first cell of a label is automatically min cell idx.
                labelMinCellIdx_.push_back(idx);
                idxStack.push(idx);
                while (!idxStack.empty()) {
                    const auto localIdx = idxStack.top();
                    idxStack.pop();

                    if (vofData.Get(localIdx, 0) > epsilon_ && labels_->GetValue(localIdx) < 0 && !isGhost(localIdx)) {
                        labels_->SetValue(localIdx, nextLabelId_);

                        const auto [i, j, k] = domainInfo_.idxToGridCoord(localIdx);

                        if (i > 0) {
                            idxStack.push(localIdx - offsetX);
                        }
                        if (i < domainInfo_.cellDims()[0] - 1) {
                            idxStack.push(localIdx + offsetX);
                        }
                        if (j > 0) {
                            idxStack.push(localIdx - offsetY);
                        }
                        if (j < domainInfo_.cellDims()[1] - 1) {
                            idxStack.push(localIdx + offsetY);
                        }
                        if (k > 0) {
                            idxStack.push(localIdx - offsetZ);
                        }
                        if (k < domainInfo_.cellDims()[2] - 1) {
                            idxStack.push(localIdx + offsetZ);
                        }
                    }
                }
                nextLabelId_++;
            }
        }
    };
} // namespace

std::tuple<std::vector<VofFlow::ComponentExtractor::BorderExtent>,
    std::vector<VofFlow::ComponentExtractor::BorderExtent>>
VofFlow::ComponentExtractor::calcBorderExtents() const {
    ZoneScoped;

    const auto& zeroExtent = domainInfo_.localZeroExtent();

    // Border cell rows within current zero extent.
    // sendExtents are at positive dimension border, receiveExtents are at negative dimension border.
    std::array<extent_t, 3> sendExtents{zeroExtent, zeroExtent, zeroExtent};
    std::array<extent_t, 3> receiveExtents{zeroExtent, zeroExtent, zeroExtent};
    sendExtents[0][0] = zeroExtent[1] - 1;
    sendExtents[1][2] = zeroExtent[3] - 1;
    sendExtents[2][4] = zeroExtent[5] - 1;
    receiveExtents[0][1] = zeroExtent[0] + 1;
    receiveExtents[1][3] = zeroExtent[2] + 1;
    receiveExtents[2][5] = zeroExtent[4] + 1;

    // Intersect all neighbor extents to find overlapping regions
    std::vector<BorderExtent> sendNeighbors;
    std::vector<BorderExtent> receiveNeighbors;

    for (int i = 0; i < 3; i++) {
        // Move extent by 1 to shift to neighbor zero extent
        auto sendNeighborExt = sendExtents[i];
        auto receiveNeighborExt = receiveExtents[i];
        sendNeighborExt[2 * i]++;
        sendNeighborExt[2 * i + 1]++;
        receiveNeighborExt[2 * i]--;
        receiveNeighborExt[2 * i + 1]--;

        for (const auto& neighbor : domainInfo_.neighbors()) {
            auto intersection = Grid::intersection(sendNeighborExt, neighbor.zeroExtent);
            if (Grid::isValidExtent(intersection)) {
                // Move intersection back to current zero extent
                intersection[2 * i]--;
                intersection[2 * i + 1]--;
                sendNeighbors.push_back({neighbor.processId, intersection});
            }
            intersection = Grid::intersection(receiveNeighborExt, neighbor.zeroExtent);
            if (Grid::isValidExtent(intersection)) {
                // Move intersection back to current zero extent
                intersection[2 * i]++;
                intersection[2 * i + 1]++;
                receiveNeighbors.push_back({neighbor.processId, intersection});
            }
        }
    }

    // Sort by process id
    std::sort(sendNeighbors.begin(), sendNeighbors.end());
    std::sort(receiveNeighbors.begin(), receiveNeighbors.end());

    return {std::move(sendNeighbors), std::move(receiveNeighbors)};
}

VofFlow::ComponentResult VofFlow::ComponentExtractor::extractComponents(
    const vtkSmartPointer<vtkDataArray>& vof) const {
    ZoneScoped;

    ExtractComponentsWorker worker(domainInfo_);

    typedef vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals> Dispatcher;
    if (!Dispatcher::Execute(vof, worker)) {
        throw std::runtime_error("Cannot dispatch array worker!");
    }

    auto labels = worker.labels_;
    auto numLabels = worker.nextLabelId_;
    auto labelMinCellIdx = std::move(worker.labelMinCellIdx_);

    if (!domainInfo_.isParallel()) {
        return {std::move(labels), numLabels};
    }

    // Idea:
    // - Generate a globally unique int64 label using process id (int32) and local label (int32).
    // - Calculate extent where neighbors are directly connected (top most cell layer at border, ignoring ghosts).
    // - Create an ordered linear list of all cell values in these regions.
    // - Send lists to neighbor, everybody sends in positive x/y/z and receives from negative x/y/z direction.
    // - Compare lists and create set of connected labels.
    // - Gather pairwise connections, create new globally unique labels and mapping from the combined ids
    // - Sort globally unique labels by global minimum cell idx for results independent of MPI process number
    // - Send mapping to everybody, and everybody updates the label names

    const int processId = mpiController_->GetLocalProcessId();
    const int numProcesses = mpiController_->GetNumberOfProcesses();

    const auto [sendNeighbors, receiveNeighbors] = calcBorderExtents();

    // Collect data to send
    std::vector<std::vector<uint64_t>> sendBuffers(sendNeighbors.size());
    std::vector<int> sendBufferSizes(sendNeighbors.size());

    for (std::size_t n = 0; n < sendNeighbors.size(); n++) {
        vtkSmartPointer<vtkIntArray> sendLabels = extractExtent(domainInfo_, sendNeighbors[n].extent, labels);
        std::vector<uint64_t> sendData(sendLabels->GetNumberOfTuples());
        for (vtkIdType i = 0; i < sendLabels->GetNumberOfTuples(); i++) {
            sendData[i] = makeCombinedLabel(processId, sendLabels->GetValue(i));
        }
        sendBufferSizes[n] = static_cast<int>(sendData.size());
        sendBuffers[n] = std::move(sendData);
    }

    // Send sizes
    std::vector<vtkMPICommunicator::Request> sendSizeRequests(sendNeighbors.size());
    for (std::size_t n = 0; n < sendNeighbors.size(); n++) {
        const int sendTo = sendNeighbors[n].neighborId;
        mpiController_->NoBlockSend(&sendBufferSizes[n], 1, sendTo, MPITags::COMPONENT_EXCHANGE_SIZE,
            sendSizeRequests[n]);
    }

    // Receive sizes
    std::vector<int> receiveBufferSizes(receiveNeighbors.size());
    std::vector<vtkMPICommunicator::Request> receiveSizeRequests(receiveNeighbors.size());
    for (std::size_t n = 0; n < receiveNeighbors.size(); n++) {
        const int receiveFrom = receiveNeighbors[n].neighborId;
        mpiController_->NoBlockReceive(&receiveBufferSizes[n], 1, receiveFrom, MPITags::COMPONENT_EXCHANGE_SIZE,
            receiveSizeRequests[n]);
    }

    // Wait
    mpiController_->WaitAll(static_cast<int>(sendSizeRequests.size()), sendSizeRequests.data());
    mpiController_->WaitAll(static_cast<int>(receiveSizeRequests.size()), receiveSizeRequests.data());

    // In ParaView 5.11.1 and on systems where long is only 32-bit (Windows), the only 64-bit type for
    // NoBlockSend/NoBlockReceive in vtkMPIController is vtkIdType. Therefore, just reinterpret the uint64_t ids
    // as vtkIdType for sending. Keep this assert, because VTK has a build option to limit vtkIdType to 32-bit.
    // ParaView 5.12 has extended the interface, so this workaround can be removed when targeting ParaView 5.12.
    static_assert(sizeof(vtkIdType) == sizeof(uint64_t), "64-bit vtkIdType is required!");

    // Send data
    std::vector<vtkMPICommunicator::Request> sendDataRequests(sendNeighbors.size());
    for (std::size_t n = 0; n < sendNeighbors.size(); n++) {
        const int sendTo = sendNeighbors[n].neighborId;
        mpiController_->NoBlockSend(reinterpret_cast<vtkIdType*>(sendBuffers[n].data()), sendBufferSizes[n], sendTo,
            MPITags::COMPONENT_EXCHANGE_DATA, sendDataRequests[n]);
    }

    // Receive data
    std::vector<std::vector<uint64_t>> receiveBuffers(receiveNeighbors.size());
    std::vector<vtkMPICommunicator::Request> receiveDataRequests(receiveNeighbors.size());
    for (std::size_t n = 0; n < receiveNeighbors.size(); n++) {
        const int receiveFrom = receiveNeighbors[n].neighborId;
        receiveBuffers[n].resize(receiveBufferSizes[n]);
        mpiController_->NoBlockReceive(reinterpret_cast<vtkIdType*>(receiveBuffers[n].data()), receiveBufferSizes[n],
            receiveFrom, MPITags::COMPONENT_EXCHANGE_DATA, receiveDataRequests[n]);
    }

    // Wait
    mpiController_->WaitAll(static_cast<int>(sendDataRequests.size()), sendDataRequests.data());
    mpiController_->WaitAll(static_cast<int>(receiveDataRequests.size()), receiveDataRequests.data());

    // Compare received data
    std::set<LabelPair> sameLabels; // TODO compare to unordered set
    for (std::size_t n = 0; n < receiveNeighbors.size(); n++) {
        vtkSmartPointer<vtkIntArray> compareLabels = extractExtent(domainInfo_, receiveNeighbors[n].extent, labels);
        for (vtkIdType i = 0; i < compareLabels->GetNumberOfTuples(); i++) {
            uint64_t receivedCombinedLabel = receiveBuffers[n][i];
            int localLabel = compareLabels->GetValue(i);
            if (getLabelIdx(receivedCombinedLabel) >= 0 && localLabel >= 0) {
                sameLabels.emplace(receivedCombinedLabel, makeCombinedLabel(processId, localLabel));
            }
        }
    }

    // Map local cell idx to global cell idx
    const auto& globalDims = domainInfo_.globalCellDims();
    const auto& gridDiff = domainInfo_.localCellDimsOffset();
    for (vtkIdType& idx : labelMinCellIdx) {
        const auto globalGridCoords = Grid::add(domainInfo_.idxToGridCoord(idx), gridDiff);
        idx = Grid::gridCoordToIdx(globalGridCoords, globalDims);
    }

    // Send all pairs to process 0 and receive local mapping
    std::vector<MappingPair> receivedMapping;
    int numLabelsGlobal = 0;

    if (processId == 0) {
        std::set<uint64_t> allLabels;
        std::vector<std::vector<vtkIdType>> allMinCellIndices;
        std::unordered_map<uint64_t, std::vector<uint64_t>> connections;

        // Process 0
        for (int l = 0; l < numLabels; l++) {
            allLabels.insert(makeCombinedLabel(0, l));
        }
        allMinCellIndices.emplace_back(labelMinCellIdx);
        for (const auto& pair : sameLabels) {
            connections[pair.label1].push_back(pair.label2);
            connections[pair.label2].push_back(pair.label1);
        }

        for (int p = 1; p < numProcesses; p++) {
            // Num labels is implicitly transferred as length of min idx list.
            std::vector<vtkIdType> labelCellIdx_p;
            ReceiveVector(mpiController_, labelCellIdx_p, p, MPITags::COMPONENT_EXCHANGE_NUM_LABELS_CELL_IDX);
            for (int l = 0; l < labelCellIdx_p.size(); l++) {
                allLabels.insert(makeCombinedLabel(p, l));
            }
            allMinCellIndices.push_back(std::move(labelCellIdx_p));
            std::vector<LabelPair> receivePairs;
            ReceiveVector(mpiController_, receivePairs, p, MPITags::COMPONENT_EXCHANGE_PAIR_LIST);
            for (const auto& pair : receivePairs) {
                connections[pair.label1].push_back(pair.label2);
                connections[pair.label2].push_back(pair.label1);
            }
        }

        // Find label groups of connected labels, store min cell idx for each group
        std::vector<std::pair<std::set<uint64_t>, vtkIdType>> labelGroups;
        while (!allLabels.empty()) {
            std::set<uint64_t> group;
            vtkIdType minCellIdx = std::numeric_limits<vtkIdType>::max();
            std::queue<uint64_t> toCheck;
            toCheck.push(*allLabels.begin());
            while (!toCheck.empty()) {
                uint64_t element = toCheck.front();
                toCheck.pop();
                group.insert(element);
                const auto& [processId, labelId] = splitCombinedLabel(element);
                minCellIdx = std::min(minCellIdx, allMinCellIndices[processId][labelId]);
                allLabels.erase(element);
                const auto it = connections.find(element);
                if (it != connections.end()) {
                    for (uint64_t conn : it->second) {
                        if (group.count(conn) == 0) {
                            toCheck.push(conn);
                        }
                    }
                    connections.erase(it); // Validation
                }
            }
            labelGroups.emplace_back(std::move(group), minCellIdx);
        }
        if (!connections.empty()) {
            // Validation
            throw std::runtime_error("Connection map not empty!");
        }

        numLabelsGlobal = static_cast<int>(labelGroups.size());

        // Sort by min cell idx for consistent naming
        std::sort(labelGroups.begin(), labelGroups.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });

        // Build mapping for each process of process local label id to global label.
        std::vector<std::vector<MappingPair>> mappings(numProcesses);
        for (std::size_t g = 0; g < labelGroups.size(); g++) {
            for (const uint64_t entry : labelGroups[g].first) {
                auto [label_process, label_id] = splitCombinedLabel(entry);
                mappings[label_process].push_back({label_id, static_cast<int>(g)});
            }
        }

        for (int p = 1; p < numProcesses; p++) {
            SendVector(mpiController_, mappings[p], p, MPITags::COMPONENT_EXCHANGE_MAPPING);
        }

        receivedMapping = std::move(mappings[0]);

    } else {
        // Num labels is implicitly transferred as length of min idx list.
        SendVector(mpiController_, labelMinCellIdx, 0, MPITags::COMPONENT_EXCHANGE_NUM_LABELS_CELL_IDX);
        std::vector<LabelPair> sendPairs(sameLabels.begin(), sameLabels.end());
        SendVector(mpiController_, sendPairs, 0, MPITags::COMPONENT_EXCHANGE_PAIR_LIST);
        ReceiveVector(mpiController_, receivedMapping, 0, MPITags::COMPONENT_EXCHANGE_MAPPING);
    }

    mpiController_->Broadcast(&numLabelsGlobal, 1, 0);

    // Create map from received data
    std::unordered_map<int, int> map;
    for (const auto& m : receivedMapping) {
        map[m.from] = m.to;
    }

    // Update local labels to global labels.
    for (vtkIdType i = 0; i < labels->GetNumberOfValues(); i++) {
        int* value = labels->GetPointer(i);
        if (*value >= 0) {
            if (domainInfo_.isInZeroExtent(i)) {
                *value = map.at(*value);
            } else {
                *value = ErrorLabels::BAD_COMPONENT_MAP;
            }
        }
    }

    return {std::move(labels), numLabelsGlobal};
}

VofFlow::ComponentResult VofFlow::ComponentExtractor::checkComponents(
    const vtkSmartPointer<vtkDataArray>& inComponents) const {
    ZoneScoped;

    auto labels = createVtkArray<vtkIntArray>(ArrayNames::LABELS, 1);
    labels->ShallowCopy(inComponents);
    labels->SetName(ArrayNames::LABELS);
    if (labels->GetNumberOfComponents() != 1 || labels->GetNumberOfTuples() != domainInfo_.numCells()) {
        throw std::runtime_error("Invalid components array!");
    }
    std::array<int, 2> range{};
    labels->GetValueRange(range.data());
    if (range[0] < ErrorLabels::EMPTY_COMPONENT) {
        throw std::runtime_error("Invalid components array with negative values!");
    }
    int numLabels = range[1] + 1;

    if (domainInfo_.isParallel()) {
        int numLabelsGlobal = 0;
        mpiController_->AllReduce(&numLabels, &numLabelsGlobal, 1, vtkCommunicator::MAX_OP);
        numLabels = numLabelsGlobal;
    }

    return {labels, numLabels};
}
