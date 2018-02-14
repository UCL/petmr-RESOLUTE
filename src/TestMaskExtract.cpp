/*
   Resolute.cpp

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

   This program tests RESOLUTE mask extraction.
 */

#include <iostream>
#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <glog/logging.h>
#include <nlohmann/json.hpp>

#include "EnvironmentInfo.h"
#include "ExtractMaskImages.hpp"

namespace fs = boost::filesystem;

int main(int argc, char **argv)
{

  //Define an input image type
  typedef itk::Image<int,4> InputImageType;

  //Define an output image type
  typedef itk::Image<int,3> OutputImageType;

  typedef ns::ExtractMaskImage<InputImageType,OutputImageType> MaskExtractorType;

  fs::path inputFile = argv[1];
  fs::path outputDir = argv[2];

  std::unique_ptr<MaskExtractorType> extractor( new MaskExtractorType(inputFile) );

  OutputImageType::Pointer frontal = OutputImageType::New();

  //extractor->GetRegion(ns::RegionLocs::FRONTAL_SINUS, frontal);
  //std::cout << frontal;

  fs::path outPath = outputDir;
  outPath /= "front-test.nii.gz";
  extractor->WriteRegionToDisk(ns::RegionLocs::FRONTAL_SINUS, outPath);

  outPath = outputDir;
  outPath /= "skull-base-test.nii.gz";
  extractor->WriteRegionToDisk(ns::RegionLocs::SKULL_BASE, outPath);

  outPath = outputDir;
  outPath /= "mastoid-test.nii.gz";
  extractor->WriteRegionToDisk(ns::RegionLocs::MASTOID, outPath);

  outPath = outputDir;
  outPath /= "nasal-test.nii.gz";
  extractor->WriteRegionToDisk(ns::RegionLocs::NASAL, outPath);

  return EXIT_SUCCESS;
}