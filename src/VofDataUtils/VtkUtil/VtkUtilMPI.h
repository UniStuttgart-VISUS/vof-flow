#pragma once

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <vtkCommunicator.h>
#include <vtkDataArray.h>
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>

namespace VofFlow {
    template<typename T>
    int AllGatherVSameLength(vtkMPIController* mpiController, const T* sendBuffer, T* recvBuffer, vtkIdType length) {
        const int numProcesses = mpiController->GetNumberOfProcesses();
        std::vector<vtkIdType> recvLengths(numProcesses, length);
        std::vector<vtkIdType> recvOffsets(numProcesses);
        for (int i = 0; i < numProcesses; ++i) {
            recvOffsets[i] = i * length;
        }
        return mpiController->AllGatherV(sendBuffer, recvBuffer, length, recvLengths.data(), recvOffsets.data());
    }

    template<typename T>
    inline int AllGatherVSameLength(vtkMPIController* mpiController, const std::vector<T>& sendBuffer,
        std::vector<T>& recvBuffer) {
        const int numProcesses = mpiController->GetNumberOfProcesses();
        if (sendBuffer.size() * numProcesses != recvBuffer.size()) {
            throw std::runtime_error("Bad buffer sizes!");
        }
        return AllGatherVSameLength(mpiController, sendBuffer.data(), recvBuffer.data(), sendBuffer.size());
    }

    template<typename T>
    int SendVector(vtkMPIController* mpiController, const std::vector<T>& data, int remoteProcessId, int tag) {
        static_assert(std::is_trivially_copyable_v<T>, "Vector type must be trivially copyable!");
        std::size_t size = data.size();
        mpiController->Send(&size, 1, remoteProcessId, tag);
        if (size > 0) {
            mpiController->Send(reinterpret_cast<const char*>(data.data()), size * sizeof(T), remoteProcessId, tag);
        }
        return 1;
    }

    template<typename T>
    int ReceiveVector(vtkMPIController* mpiController, std::vector<T>& data, int remoteProcessId, int tag) {
        static_assert(std::is_trivially_copyable_v<T>, "Vector type must be trivially copyable!");
        std::size_t size = 0;
        mpiController->Receive(&size, 1, remoteProcessId, tag);
        if (size > 0) {
            data.resize(size);
            mpiController->Receive(reinterpret_cast<char*>(data.data()), size * sizeof(T), remoteProcessId, tag);
        } else {
            data.clear();
        }
        return 1;
    }

    template<typename T>
    T ReduceSum(vtkMPIController* mpiController, T value, int destProcessId = 0) {
        T globalValue{};
        mpiController->Reduce(&value, &globalValue, 1, vtkCommunicator::StandardOperations::SUM_OP, destProcessId);
        return globalValue;
    }

    inline int ReduceArray(vtkMPIController* mpiController, vtkDataArray* buffer, int operation, int destProcessId) {
        // MPI reduce requires an independent receive buffer. In addition, all metadata in the `buffer` array should
        // not change. We use DeepCopy here to copy full metadata. The MPI controller will allocate the receive-buffer
        // anyway, so the only overhead is the memory copy. Alternative would be to either manually copy all metadata
        // which requires knowing all possible metadata or to copy all values after the reduce to the original `buffer`
        // array which would require an array dispatch worker to avoid using a loop with the generic double typed
        // Get/Set interface for individual values.
        // As we only require output for the process `destProcessId` we can use a dummy receive-buffer without correct
        // metadata on all other processes. For `destProcessId`, we can copy the array before reduce and use the
        // original array directly as receive buffer.
        int status = 0;
        if (mpiController->GetLocalProcessId() == destProcessId) {
            vtkSmartPointer<vtkDataArray> sendBuffer = buffer->NewInstance();
            sendBuffer->DeepCopy(buffer);
            status = mpiController->Reduce(sendBuffer, buffer, operation, destProcessId);
        } else {
            vtkSmartPointer<vtkDataArray> recvBuffer = buffer->NewInstance();
            recvBuffer->Initialize();
            status = mpiController->Reduce(buffer, recvBuffer, operation, destProcessId);
        }
        return status;
    }

} // namespace VofFlow
