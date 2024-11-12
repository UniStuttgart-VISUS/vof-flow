#include "VtkUtilArray.h"

#include <vtkAOSDataArrayTemplate.h>
#include <vtkArrayDispatch.h>

namespace {
    struct ShallowCopyWorker {
        vtkSmartPointer<vtkDataArray> copy_;

        template<typename ValueType>
        void operator()(vtkAOSDataArrayTemplate<ValueType>* dataArray) {
            copy_ = vtkSmartPointer<vtkAOSDataArrayTemplate<ValueType>>::New();
            copy_->ShallowCopy(dataArray);
        }
    };
} // namespace

vtkSmartPointer<vtkDataArray> VofFlow::shallowCopy(const vtkSmartPointer<vtkDataArray>& data) {
    ShallowCopyWorker worker;
    typedef vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals> Dispatcher;
    if (!Dispatcher::Execute(data, worker)) {
        throw std::runtime_error("Cannot dispatch array worker!");
    }
    return worker.copy_;
}
