/*
   ExtractMaskImages.hpp

   Author:      Benjamin A. Thomas

   Copyright 2018 Institute of Nuclear Medicine, University College London.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
 */

#pragma once

#ifndef _EXTRACTMASKIMAGES_H_
#define _EXTRACTMASKIMAGES_H_

#include <boost/filesystem.hpp>
#include <glog/logging.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageDuplicator.h>
#include <itkExtractImageFilter.h>

namespace ns {

enum class RegionLocs {
FRONTAL_SINUS = 0,
SKULL_BASE = 1,
MASTOID = 2,
NASAL = 3
};

template <class TImage, class TOutImage>
class ExtractMaskImage
{

public:

  typedef TImage ImageType;
  typedef TOutImage OutputImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typedef itk::ExtractImageFilter< ImageType, OutputImageType > ExtractFilterType;

  ExtractMaskImage(boost::filesystem::path inputFilePath);

  void GetRegion(const RegionLocs regID, typename OutputImageType::Pointer &dst );
  void WriteRegionToDisk(const RegionLocs regID, const boost::filesystem::path outFileName );

protected:

   typename ImageType::Pointer _inputImage;


};

template <typename TImage, class TOutImage>
ExtractMaskImage<TImage, TOutImage>::ExtractMaskImage(boost::filesystem::path inputFilePath){

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputFilePath.string().c_str());

  try {
    reader->Update();
    //Duplicate contents of reader into _inputImage.
    typedef itk::ImageDuplicator<ImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(reader->GetOutput());
    duplicator->Update();
    _inputImage = duplicator->GetOutput();
  } catch (itk::ExceptionObject &ex) {
    LOG(ERROR) << ex;
    LOG(ERROR) << "Unable to get image from file: " << inputFilePath;
    throw false;
  }

  LOG(INFO) << "Mask dimensions: " << _inputImage->GetLargestPossibleRegion().GetSize();


}

template <typename TImage, class TOutImage>
void ExtractMaskImage<TImage, TOutImage>::GetRegion(RegionLocs regID, typename TOutImage::Pointer &dst){

  //Extract filter used to extract nd-1 volume from an nd file.
  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  extractFilter->SetInput(_inputImage);
  extractFilter->SetDirectionCollapseToIdentity(); 

  //Get mask image size.
  typename ImageType::SizeType imageSize =
    _inputImage->GetLargestPossibleRegion().GetSize();

  typename ImageType::IndexType desiredStart;
  desiredStart.Fill(0);
  typename ImageType::SizeType desiredSize = imageSize;

  //Starts reading from 4D volume at index (0,0,0,i) through to
  //(maxX, maxY, maxZ,0).
  desiredStart[3] = (int)regID;
  desiredSize[3] = 0;

  typename ImageType::RegionType region(desiredStart, desiredSize);

  extractFilter->SetExtractionRegion( region );

  try {
    extractFilter->Update();
    
    //Duplicate contents of reader into dst.
    typedef itk::ImageDuplicator<OutputImageType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(extractFilter->GetOutput());
    duplicator->Update();
    dst = duplicator->GetOutput();

  } catch (itk::ExceptionObject &ex) {
    LOG(ERROR) << ex;
    LOG(ERROR) << "Unable to extract region!";
    throw false;
  }

}

template <typename TImage, class TOutImage>
void ExtractMaskImage<TImage, TOutImage>::WriteRegionToDisk(
  const RegionLocs regID, 
  const boost::filesystem::path outFileName ){

  typename OutputImageType::Pointer outImage;

  try {
    GetRegion( regID, outImage);
  }
  catch (bool) {
    LOG(ERROR) << "Extraction failed!";
    throw false;
  }

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outFileName.string().c_str());

  try {
    writer->SetInput(outImage);
    writer->Update();
  } catch (itk::ExceptionObject &ex) {
    LOG(ERROR) << ex;
    LOG(ERROR) << "Unable to extract region!";
    throw false;
  }
  
}

} //end namespace ns

#endif 
