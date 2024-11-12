#include "VtkUtilWriter.h"

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLRectilinearGridWriter.h>

void VofFlow::writeData(vtkPolyData* data, const std::string& path, bool compressed) {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetInputData(data);
    writer->SetFileName(path.c_str());
    writer->SetDataModeToAppended();
    writer->SetEncodeAppendedData(false);
    if (compressed) {
        writer->SetCompressorTypeToZLib();
        writer->SetCompressionLevel(1);
    } else {
        writer->SetCompressorTypeToNone();
    }
    writer->SetHeaderTypeToUInt64();
    writer->Write();
}

void VofFlow::writeData(vtkRectilinearGrid* data, const std::string& path, bool compressed) {
    vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
    writer->SetInputData(data);
    writer->SetFileName(path.c_str());
    writer->SetDataModeToAppended();
    writer->SetEncodeAppendedData(false);
    if (compressed) {
        writer->SetCompressorTypeToZLib();
        writer->SetCompressionLevel(1);
    } else {
        writer->SetCompressorTypeToNone();
    }
    writer->SetHeaderTypeToUInt64();
    writer->Write();
}

void VofFlow::writeData(vtkImageData* data, const std::string& path, bool compressed) {
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetInputData(data);
    writer->SetFileName(path.c_str());
    writer->SetDataModeToAppended();
    writer->SetEncodeAppendedData(false);
    if (compressed) {
        writer->SetCompressorTypeToZLib();
        writer->SetCompressionLevel(1);
    } else {
        writer->SetCompressorTypeToNone();
    }
    writer->SetHeaderTypeToUInt64();
    writer->Write();
}
