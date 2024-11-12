#pragma once

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkType.h>

#include "../Math/Vector.h"

namespace VofFlow {
    vtkSmartPointer<vtkDataArray> shallowCopy(const vtkSmartPointer<vtkDataArray>& data);

    template<typename ArrayT>
    inline vtkSmartPointer<ArrayT> createVtkArray(const char* name, int numComponents = 1, vtkIdType numTuples = 0) {
        vtkSmartPointer<ArrayT> a = vtkSmartPointer<ArrayT>::New();
        a->SetName(name);
        a->SetNumberOfComponents(numComponents);
        if (numTuples > 0) {
            a->SetNumberOfTuples(numTuples);
        }
        return a;
    }

    template<typename ArrayT>
    inline vtkSmartPointer<ArrayT> createVtkArray(const char* name, int numComponents, vtkIdType numTuples,
        double fill) {
        auto a = createVtkArray<ArrayT>(name, numComponents, numTuples);
        a->Fill(fill);
        return a;
    }

    template<typename ArrayT>
    inline vtkSmartPointer<ArrayT> createVtkArray(const char* name,
        const std::vector<typename ArrayT::ValueType>& data) {
        auto a = createVtkArray<ArrayT>(name, 1, static_cast<vtkIdType>(data.size()));
        std::memcpy(a->GetVoidPointer(0), data.data(), sizeof(typename ArrayT::ValueType) * data.size());
        return a;
    }

    inline vtkSmartPointer<vtkFloatArray> createVtkArray(const char* name, const std::vector<vec3>& data) {
        auto a = createVtkArray<vtkFloatArray>(name, 3, static_cast<vtkIdType>(data.size()));
        std::memcpy(a->GetVoidPointer(0), data.data(), sizeof(vec3) * data.size());
        return a;
    }

    inline vtkSmartPointer<vtkStringArray> createVtkStringArrayFromString(const char* name, const std::string& str) {
        vtkSmartPointer<vtkStringArray> a = vtkSmartPointer<vtkStringArray>::New();
        a->SetName(name);
        a->SetNumberOfValues(1);
        a->SetValue(0, str);
        return a;
    }

    template<typename ArrayT>
    void appendArray(ArrayT* a, ArrayT* b) {
        if (a->GetNumberOfComponents() != b->GetNumberOfComponents()) {
            throw std::runtime_error("Bad array format!");
        }
        auto* dest = a->WritePointer(a->GetNumberOfValues(), b->GetNumberOfValues());
        const auto* src = b->GetPointer(0);
        std::memcpy(dest, src, b->GetNumberOfValues() * sizeof(typename ArrayT::ValueType));
    }

    template<typename ArrayT>
    void appendArray(vtkSmartPointer<ArrayT>& a, vtkSmartPointer<ArrayT>& b) {
        appendArray(a.Get(), b.Get());
    }

    template<typename ArrayT>
    inline typename ArrayT::ValueType* getTuplePointer(ArrayT* array, vtkIdType tupleIdx) {
        return array->GetPointer(tupleIdx * array->GetNumberOfComponents());
    }

    template<typename ArrayT>
    inline typename ArrayT::ValueType* getTuplePointer(const vtkSmartPointer<ArrayT>& array, vtkIdType tupleIdx) {
        return getTuplePointer(array.Get(), tupleIdx);
    }

    template<typename T>
    std::vector<T> getArrayValues(vtkDataArray* data) {
        const vtkIdType size = data->GetNumberOfTuples();
        std::vector<T> v(size);
        for (vtkIdType i = 0; i < size; i++) {
            v[i] = static_cast<T>(data->GetComponent(i, 0));
        }
        return v;
    }

    inline void appendArrayName(vtkDataArray* data, const std::string& nameExtra) {
        data->SetName((std::string(data->GetName()) + nameExtra).c_str());
    }

    inline vec3& getVec3TuplePointer(vtkFloatArray* array, vtkIdType tupleIdx) {
        // assume NumberOfComponents == 3
        return *reinterpret_cast<vec3*>(getTuplePointer(array, tupleIdx));
    }

} // namespace VofFlow
