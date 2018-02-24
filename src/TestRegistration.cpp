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

#include "antsRegistrationTemplateHeader.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

namespace fs = boost::filesystem;
using namespace ants;

int main(int argc, char **argv)
{

std::string argList = "--verbose 1 --dimensionality 3 \
--float 0 --collapse-output-transforms 1 \
--output [reg-out/fullSyN,reg-out/fullSyNWarped.nii.gz,reg-out/fullSyNInverseWarped.nii.gz] \
--interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
--initial-moving-transform [/Users/bathomas/Documents/ATLASES/mni_icbm152_nlin_sym_09a_nifti/mni_icbm152_nlin_sym_09a/mni_icbm152_t1_tal_nlin_sym_09a.nii,/Users/bathomas/Documents/BUILDS/PETMR-RESOLUTE/ute2-test.nii.gz,1] \
--transform Rigid[0.1] \
--metric MI[/Users/bathomas/Documents/ATLASES/mni_icbm152_nlin_sym_09a_nifti/mni_icbm152_nlin_sym_09a/mni_icbm152_t1_tal_nlin_sym_09a.nii,/Users/bathomas/Documents/BUILDS/PETMR-RESOLUTE/ute2-test.nii.gz,1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
--transform Affine[0.1] \
--metric MI[/Users/bathomas/Documents/ATLASES/mni_icbm152_nlin_sym_09a_nifti/mni_icbm152_nlin_sym_09a/mni_icbm152_t1_tal_nlin_sym_09a.nii,/Users/bathomas/Documents/BUILDS/PETMR-RESOLUTE/ute2-test.nii.gz,1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
--transform BSplineSyN[0.1,26,0,3] \
--metric CC[/Users/bathomas/Documents/ATLASES/mni_icbm152_nlin_sym_09a_nifti/mni_icbm152_nlin_sym_09a/mni_icbm152_t1_tal_nlin_sym_09a.nii,/Users/bathomas/Documents/BUILDS/PETMR-RESOLUTE/ute2-test.nii.gz,1,4] \
--convergence [10x7x5x5,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox";

std::vector<std::string> args;

//const boost::regex ex("--");
//boost::split(args, argList, ex, boost::token_compress_on);

boost::split_regex( args, argList, boost::regex( " " ) ) ;

int x=1;
for (auto a : args){
  LOG(INFO) << x << "\t\t" << a;
  x++;
}

/*
args.push_back("verbose 1");
args.push_back("dimensionality 3");
args.push_back("float 0");
args.push_back("collapse-output-transforms 1");
*/

antsRegistration( args, &std::cout);


  /*
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
  extractor->WriteRegionToDisk(ns::RegionLocs::NASAL, outPath);*/

  return EXIT_SUCCESS;
}