#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <vtkCommunicator.h>
#include <vtkMPIController.h>
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
} // namespace VofFlow
