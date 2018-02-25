/*
   ANTsReg.hpp

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

#ifndef _ANTSREG_HPP_
#define _ANTSREG_HPP_ 

#include <nlohmann/json.hpp>
#include <boost/filesystem.hpp>
#include <glog/logging.h>

#include <string>

#include <itkImageFileReader.h>

#include <antsRegistrationTemplateHeader.h>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>


namespace reg {

template <class TImage>
class ANTsReg
{

public:
  typedef TImage ImageType;

  //DicomReader();
  ANTsReg();

  void SetParams(const nlohmann::json &params) { _regParams = params; };
  void SetOutputDirectory(const boost::filesystem::path outDir);
  void SetOutputPrefix(const std::string s){ _prefix = s; };
  void SetReferenceFileName(const boost::filesystem::path refFileName);
  void SetFloatingFileName(const boost::filesystem::path floatFileName);

  void Update();

  typename TImage::ConstPointer GetOutputImage();
  typename TImage::ConstPointer GetOutputInverseImage();

protected:

  void InsertParam(const std::string &key, const std::string &info);

  std::string _prefix;
  boost::filesystem::path _outDir;
  boost::filesystem::path _refFileName;
  boost::filesystem::path _floatFileName;

  nlohmann::json _regParams;
  std::string _argList;

};

//Constructor
template <typename TImage>
ANTsReg<TImage>::ANTsReg()
{

  DLOG(INFO) << "Initialised ANTsReg.";

  _argList = "--verbose 1 --dimensionality 3 \
  --float 0 --collapse-output-transforms 1 \
  --output [<%%PREFIX%%>,<%%WARPEDIMG%%>,<%%INVWARPEDIMG%%>] \
  --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
  --initial-moving-transform [<%%REF%%>,<%%FLOAT%%>,1] \
  --transform Rigid[0.1] \
  --metric MI[<%%REF%%>,<%%FLOAT%%>,1,32,Regular,0.25] \
  --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
  --transform Affine[0.1] \
  --metric MI[<%%REF%%>,<%%FLOAT%%>,1,32,Regular,0.25] \
  --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
  --transform BSplineSyN[0.1,26,0,3] \
  --metric CC[<%%REF%%>,<%%FLOAT%%>,1,4] \
  --convergence [10x7x5x5,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox";

}

template <typename TImage>
typename TImage::ConstPointer ANTsReg<TImage>::GetOutputImage(){

  boost::filesystem::path tImagePath = _outDir;
  tImagePath /= _prefix;
  tImagePath += "Warped.nii.gz";

  typedef typename itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  try {
    reader->SetFileName( tImagePath.string() );
    reader->Update();
    LOG(INFO) << "Read " << tImagePath << " successfully.";
  }
  catch (itk::ExceptionObject &ex)
  {
    //std::cout << ex << std::endl;
    LOG(ERROR) << ex;
    LOG(ERROR) << "Unable to read image " << tImagePath;
    throw false;
  }

  return static_cast< const TImage * >( reader->GetOutput() );
}

template <typename TImage>
typename TImage::ConstPointer ANTsReg<TImage>::GetOutputInverseImage(){

  boost::filesystem::path iImagePath = _outDir;
  iImagePath /= _prefix;
  iImagePath += "InverseWarped.nii.gz";

  typedef typename itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  try {
    reader->SetFileName( iImagePath.string() );
    reader->Update();
    LOG(INFO) << "Read " << iImagePath << " successfully.";
  }
  catch (itk::ExceptionObject &ex)
  {
    //std::cout << ex << std::endl;
    LOG(ERROR) << ex;
    LOG(ERROR) << "Unable to read image " << iImagePath;
    throw false;
  }

  return static_cast< const TImage * >( reader->GetOutput() );
}

template <typename TImage>
void ANTsReg<TImage>::SetOutputDirectory(const boost::filesystem::path outDir){

  if ( boost::filesystem::exists(outDir) ){
    if (!boost::filesystem::is_directory(outDir)){
      LOG(ERROR) << "Output directory: " << outDir << "does not appear to be a directory!";
      throw(false);
    }
  } else {
    try {
      boost::filesystem::create_directories(outDir);
    } catch (const boost::filesystem::filesystem_error &e){
      LOG(ERROR) << "Cannot create output folder : " << outDir;
      throw(false);
    }
  }

  _outDir = outDir;

}

template <typename TImage>
void ANTsReg<TImage>::SetReferenceFileName(const boost::filesystem::path refFileName){

  try {
    if (!boost::filesystem::is_regular_file(refFileName)){
      LOG(ERROR) << "Cannot find reference image file : " << refFileName;
      throw(false);      
    }
  } catch (const boost::filesystem::filesystem_error &e){
      LOG(ERROR) << "File system error reading : " << refFileName;
      throw(false);
    }
  
  _refFileName = refFileName;
  
}

template <typename TImage>
void ANTsReg<TImage>::SetFloatingFileName(const boost::filesystem::path floatFileName){

  try {
    if (!boost::filesystem::is_regular_file(floatFileName)){
      LOG(ERROR) << "Cannot find floating image file : " << floatFileName;
      throw(false);      
    }
  } catch (const boost::filesystem::filesystem_error &e){
      LOG(ERROR) << "File system error reading : " << floatFileName;
      throw(false);
    }
  
  _floatFileName = floatFileName;
  
}

template <typename TImage>
void ANTsReg<TImage>::InsertParam(const std::string &key, const std::string &info){


  std::string target = "<%%" + key + "%%>";
  std::string::size_type n = _argList.find(target);

  if (n == std::string::npos){
    LOG(WARNING) << "Replacement key: " << target << " not found!";
  }
  //_argList.replace(n,target.length(),info);
  boost::replace_all(_argList, target, info);
}

template <typename TImage>
void ANTsReg<TImage>::Update(){

/*std::string argList = "--verbose 1 --dimensionality 3 \
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
--convergence [10x7x5x5,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox";*/

  boost::filesystem::path fullPrefix = _outDir;
  fullPrefix /= _prefix;

  boost::filesystem::path tImagePath = fullPrefix;
  tImagePath += "Warped.nii.gz";

  boost::filesystem::path iImagePath = fullPrefix;
  iImagePath += "InverseWarped.nii.gz";

  InsertParam("FLOAT", _floatFileName.string());
  InsertParam("REF", _refFileName.string());
  InsertParam("PREFIX", fullPrefix.string());
  InsertParam("WARPEDIMG", tImagePath.string());
  InsertParam("INVWARPEDIMG", iImagePath.string());

  std::vector<std::string> args;

//const boost::regex ex("--");
//boost::split(args, argList, ex, boost::token_compress_on);

  boost::split_regex( args, _argList, boost::regex( " " ) ) ;

  std::vector<std::string> finalArgs;

  int x=1;
  for (auto a : args){
    if (a != "") {
      finalArgs.push_back(a);
      LOG(INFO) << x << "\t\t" << a;
      x++;
    }
  }

  LOG(INFO) << "Starting ANTs registration. This will take a while...";
    ants::antsRegistration( finalArgs, &std::cout);
  LOG(INFO) << "Registration complete!";

  
}

} //namespace reg


#endif