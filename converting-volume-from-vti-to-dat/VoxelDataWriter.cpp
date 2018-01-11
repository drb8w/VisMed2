#include "VoxelDataWriter.h"

#include <vtkImageResample.h>
#include <vtkAlgorithmOutput.h>

#include <algorithm>

bool VoxelDataWriter::WriteDatFile(vtkSmartPointer<vtkImageData> imageData, const std::string &voxelDataFileName, bool invertValue)
{
	// write  12bit, row major dat file for reading with visualization app
	// every entry consists of an unsigned short, so 16 bit

	FILE* pFile;
	pFile = fopen(voxelDataFileName.c_str(), "wb");

	if (pFile != nullptr)
	{
		// --------------------------------------
		// https://public.kitware.com/pipermail/vtkusers/2009-November/054596.html

		int* dims = imageData->GetDimensions();

		unsigned short width = (unsigned short)dims[0];
		unsigned short height = (unsigned short)dims[1];
		unsigned short depth = (unsigned short)dims[2];

		fwrite(&width, sizeof(unsigned short), 1, pFile);
		fwrite(&height, sizeof(unsigned short), 1, pFile);
		fwrite(&depth, sizeof(unsigned short), 1, pFile);

		// <DataArray type="Int16" Name="Scalars_" format="appended" RangeMin="-1024" RangeMax="3010" offset="0"/>

		int type = imageData->GetScalarType();
		double min = imageData->GetScalarTypeMin();
		double max = imageData->GetScalarTypeMax();

		double range[2];
		imageData->GetScalarRange(range);

		for (int z = 0; z < dims[2]; z++)
		{
			for (int y = 0; y < dims[1]; y++)
			{
				for (int x = 0; x < dims[0]; x++)
				{
					// zero is the component, add another loop if you have more
					// than one component
					double intensity = (imageData->GetScalarComponentAsDouble(x, y, z, 0) - range[0]) / (range[1] - range[0]);
					// do something with intensity
					if (invertValue)
						intensity = 1.0 - intensity;
					unsigned short us = (unsigned short)(intensity * (pow(2, 12) - 1));
					fwrite(&us, sizeof(unsigned short), 1, pFile);
				}
			}
		}

		// --------------------------------------

		fclose(pFile);

		return true;
	}

	return false;

	// =========================================================================================

	//MV_RGBA color = pVoxelModel->palette[pData[idx_z + idx_y * stride_y + idx_x * stride_x]];

	//// https://stackoverflow.com/questions/687261/converting-rgb-to-grayscale-intensity
	//// https://en.wikipedia.org/wiki/Grayscale
	//// RGB scaling: 0.2126 * R + 0.7152 * G + 0.0722 * B
	////// not premultiplied alpha
	////return (color.r * 0.2126f + color.g * 0.7152f + color.b * 0.0722f) * color.a / pow(255.0f, 2.0f);
	//// premultiplied alpha
	//float intensity = (color.r * 0.2126f + color.g * 0.7152f + color.b * 0.0722f) / 255.0f;
	//if (invertValue)
	//	intensity = 1.0f - intensity;
	//unsigned short us = (unsigned short)(intensity * (pow(2, 12) - 1));
	//fwrite(&us, sizeof(unsigned short), 1, pFile);
}

bool VoxelDataWriter::WriteDatFile(vtkSmartPointer<vtkXMLImageDataReader> reader, const std::string &voxelDataFileName, bool invertValue, unsigned short resolution_max)//, bool blockScale)
{
	// write  12bit, row major dat file for reading with visualization app
	// every entry consists of an unsigned short, so 16 bit

	// --------------------------------------
	// read image
	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	imageData->ShallowCopy(reader->GetOutput());

	int* dims = imageData->GetDimensions();
	// calculate scaling factors
	double factors[3];
	
	//if (blockScale)
	//{
	//	// find maximum
	//	double size_max = dims[0];
	//	if (size_max < dims[1])
	//		size_max = dims[1];
	//	if (size_max < dims[2])
	//		size_max = dims[2];

	//	// calculate scaling factors
	//	factors[0] = size_max / dims[0];
	//	factors[1] = size_max / dims[1];
	//	factors[2] = size_max / dims[2];
	//}
	//else 
	//{
		imageData->ComputeBounds();

		// https://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/VolumeRendering/Cxx/FixedPointVolumeRayCastMapperCT.cxx
		// (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
		double *bounds = imageData->GetBounds();
		//int* dims = imageData->GetDimensions();

		// resample to get squared pixels
		double size_vox[3];
		size_vox[0] = (bounds[1] - bounds[0]) / dims[0];
		size_vox[1] = (bounds[3] - bounds[2]) / dims[1];
		size_vox[2] = (bounds[5] - bounds[4]) / dims[2];

		// find maximum
		double size_min = size_vox[0];
		if (size_min > size_vox[1])
			size_min = size_vox[1];
		if (size_min > size_vox[2])
			size_min = size_vox[2];

		// calculate scaling factors
		factors[0] = size_vox[0] / size_min;
		factors[1] = size_vox[1] / size_min;
		factors[2] = size_vox[2] / size_min;

		// croping to maximal resolution
		double dim_max = std::max(std::max(factors[0] * dims[0], factors[1] * dims[1]), factors[2] * dims[2]);
	
		if (dim_max > resolution_max) {
			double scale_max = resolution_max / dim_max;
			factors[0] *= scale_max;
			factors[1] *= scale_max;
			factors[2] *= scale_max;
		}
	//}

	vtkImageResample *resample = vtkImageResample::New();
	resample->SetInputConnection(reader->GetOutputPort());
	resample->SetAxisMagnificationFactor(0, factors[0]);
	resample->SetAxisMagnificationFactor(1, factors[1]);
	resample->SetAxisMagnificationFactor(2, factors[2]);

	resample->Update();

	// --------------------------------------
	vtkSmartPointer<vtkImageData> imageData2 = vtkSmartPointer<vtkImageData>::New();
	imageData2->ShallowCopy(resample->GetOutput());
	// --------------------------------------

	return WriteDatFile(imageData2, voxelDataFileName, invertValue);
}

bool VoxelDataWriter::WriteDatFile(vtkSmartPointer<vtkXMLImageDataReader> reader, double *scaleFactors, const std::string &voxelDataFileName, bool invertValue)
{
	// write  12bit, row major dat file for reading with visualization app
	// every entry consists of an unsigned short, so 16 bit

	// --------------------------------------
	// read image
	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	imageData->ShallowCopy(reader->GetOutput());

	vtkImageResample *resample = vtkImageResample::New();
	resample->SetInputConnection(reader->GetOutputPort());
	resample->SetAxisMagnificationFactor(0, scaleFactors[0]);
	resample->SetAxisMagnificationFactor(1, scaleFactors[1]);
	resample->SetAxisMagnificationFactor(2, scaleFactors[2]);

	resample->Update();

	// --------------------------------------
	vtkSmartPointer<vtkImageData> imageData2 = vtkSmartPointer<vtkImageData>::New();
	imageData2->ShallowCopy(resample->GetOutput());
	// --------------------------------------

	return WriteDatFile(imageData2, voxelDataFileName, invertValue);

}

