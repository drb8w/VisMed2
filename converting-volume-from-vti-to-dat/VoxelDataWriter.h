#pragma once

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>

class VoxelDataWriter {

public:
	static bool WriteDatFile(vtkSmartPointer<vtkImageData> imageData, const std::string &voxelDataFileName, bool invertValue);
	
	//static bool WriteDatFile(vtkSmartPointer<vtkXMLImageDataReader> reader, const std::string &voxelDataFileName, bool invertValue);
	//static bool WriteDatFile(vtkSmartPointer<vtkXMLImageDataReader> reader, const std::string &voxelDataFileName, bool invertValue, bool blockScale = false);
	static bool WriteDatFile(vtkSmartPointer<vtkXMLImageDataReader> reader, const std::string &voxelDataFileName, bool invertValue, unsigned short resolution_max);// , bool blockScale = false);

	static bool WriteDatFile(vtkSmartPointer<vtkXMLImageDataReader> reader, double *scaleFactors, const std::string &voxelDataFileName, bool invertValue);

};
