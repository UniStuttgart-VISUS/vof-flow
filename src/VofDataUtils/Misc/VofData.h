#pragma once

#include <vtkDataArray.h>
#include <vtkSmartPointer.h>

namespace VofFlow {
    /*
     * In two phases/species data there is only a single field usually called `vof` (`funs*.hdf`). In three
     * phases/species data there are two fields called `vof` (`funs*.hdf`) and `f3` (`fun3*.hdf`). Here, the PLIC
     * surface of the `f3` field is based on the gradient, while the PLIC surface of the `vof` field uses additional
     * normals in the dataset `n_c_3ph` (`n3ph*.hdf`). The `f3` field in the three phases/species case and the `vof`
     * field in the two phases/species case are handled the same way, despite the name being different.
     * Therefore, we use the terms `vof1st` and `vof2nd` to simplify the different cases. So `vof1st` is always the
     * field using the gradient, while `vof2nd` (if present) is the grid using the additional normals.
     */
    struct VofData {
        vtkSmartPointer<vtkDataArray> vof1st;
        vtkSmartPointer<vtkDataArray> vof2nd;
        vtkSmartPointer<vtkDataArray> vof2ndNorm;
    };
} // namespace VofFlow
