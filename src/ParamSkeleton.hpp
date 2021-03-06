/*
   ParamSkeleton.hpp

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

#ifndef _PARAMSKELETON_HPP_
#define _PARAMSKELETON_HPP_ 

#include <nlohmann/json.hpp>
#include <glog/logging.h>

namespace ns {

  struct params {

    std::string version;

    std::string destFileType;
    boost::filesystem::path destDir;
    std::string destExportMethod;

    boost::filesystem::path logDir;

    std::string MRACSeriesName;
    std::string UTE1SeriesName;
    std::string UTE1TE;
    std::string UTE2SeriesName;
    std::string UTE2TE;

    std::string regName;
    boost::filesystem::path regTemplatePath;
    std::string regArgs;
  };

  void to_json(nlohmann::json &j, const params &p){

    j=nlohmann::json{ 
        {"version", p.version},
        {"destFileType", p.destFileType},
        {"destDir", p.destDir.string()},
        {"destExportMethod", p.destExportMethod},

        {"logDir", p.logDir.string()},

        {"MRACSeriesName", p.MRACSeriesName},
        {"UTE1SeriesName", p.UTE1SeriesName},
        {"UTE1TE", p.UTE1TE},
        {"UTE2SeriesName", p.UTE2SeriesName},
        {"UTE2TE", p.UTE2TE},

        {"regName", p.regName},
        {"regTemplatePath", p.regTemplatePath.string()},
        {"regArgs", p.regArgs}
    };
  }

  void from_json(const nlohmann::json &j, params &p){
    //LOG(INFO) << j.dump(4);
    p.version = j.at("version").get<std::string>();
    p.destFileType = j.at("destFileType").get<std::string>();
    p.destDir = j.at("destDir").get<std::string>();
    p.destExportMethod = j.at("destExportMethod").get<std::string>();
    p.logDir = j.at("logDir").get<std::string>();

    p.MRACSeriesName = j.at("MRACSeriesName").get<std::string>();
    p.UTE1SeriesName = j.at("UTE1SeriesName").get<std::string>();
    p.UTE1TE = j.at("UTE1TE").get<std::string>();   
    p.UTE2SeriesName = j.at("UTE2SeriesName").get<std::string>();
    p.UTE2TE = j.at("UTE2TE").get<std::string>();    

    p.regName = j.at("regName").get<std::string>();
    p.regTemplatePath = j.at("regTemplatePath").get<std::string>();
    p.regArgs = j.at("regArgs").get<std::string>();

  }

  const params skeleton = {
    VERSION_NO,
    ".nii.gz",
    ".",
    "FILE",

    "./logs",

    "Head_MRAC_PET_UTE_UMAP",
    "Head_MRAC_PET_UTE",
    "0.07",
    "Head_MRAC_PET_UTE",
    "2.46",

    "ANTS",
    "",
    "3 -m CC[<%%REF%%>,<%%FLOAT%%>,1,4] -i 10x5x2 -o <%%PREFIX%%> -t SyN[0.5] -r Gauss[3,0] -G"
    //"--verbose 0 --dimensionality 3 --float 1 --collapse-output-transforms 1 --output [<%%PREFIX%%>,<%%WARPEDIMG%%>,<%%INVWARPEDIMG%%>] --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] --initial-moving-transform [<%%REF%%>,<%%FLOAT%%>,1] --transform Affine[0.1] --metric MI[<%%REF%%>,<%%FLOAT%%>,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[0.5,3,0] --metric CC[<%%REF%%>,<%%FLOAT%%>,1,4] --convergence [10x5x2,1e-6,10] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0mm",
  };

bool ValidateJSON(const nlohmann::json j){

  ns::params pr;

  //If there is an error in the JSON file.
  try {
    pr = j;
  }
  catch (nlohmann::json::exception& e){
    LOG(ERROR) << e.what();
    throw false;
  }

  //Try to make the log directory if it doesn't exist.
  if (!boost::filesystem::exists(pr.logDir)){
    LOG(WARNING) << "Logging destination does not exist";
    LOG(INFO) << "Attempting to create logging directory...";
    try {
      boost::filesystem::create_directories(pr.logDir);
    } catch (const boost::filesystem::filesystem_error &e) {
      LOG(ERROR) << " cannot create folder!";
      throw false;
    }
  }

  //If the log destination given in the JSON file is a file, then throw
  if (boost::filesystem::is_regular_file(pr.logDir)){
      LOG(ERROR) << "LogPath is a file. This should be a directory!";
      throw false;    
  }

  //Check if registration template exists
  if (!(boost::filesystem::exists(pr.regTemplatePath.string())))  {
    LOG(ERROR) << "Registration template " << pr.regTemplatePath.string() << " does not exist!";
    throw false;
  }

  //Check if registration template is a file
  if (!boost::filesystem::is_regular_file(pr.regTemplatePath)){
    LOG(ERROR) << "Registration template " << pr.regTemplatePath << " does not appear to be a file!";
    throw false;
  }

  /*
  //Check if registration mask exists
  if (!boost::filesystem::exists(pr.regMaskPath)){
    LOG(ERROR) << "Registration mask " << pr.regMaskPath << " does not exist!";
    throw false;
  }

  //Check if registration mask is a file
  if (!boost::filesystem::is_regular_file(pr.regMaskPath)){
    LOG(ERROR) << "Registration mask " << pr.regMaskPath << " does not appear to be a file!";
    throw false;
  }*/

  return true;
}

void WriteJSONSkeleton(const boost::filesystem::path &outFile){
  //Write an initial JSON parameter file to outFile.

  DLOG(ERROR) << "outFile = " << outFile;

  std::string dst = outFile.string();
  boost::algorithm::replace_all(dst, "\"", "");

  if (boost::filesystem::exists(dst)){
    LOG(ERROR) << "JSON config. file already exists!";
    LOG(ERROR) << "Refusing to overwrite!";
    throw false;
  }

  std::ofstream file(dst);

  if (!file){
    LOG(ERROR) << "Cannot open output file for writing!" << std::endl; 
    throw false;
  }

  nlohmann::json jsonSkeleton = ns::skeleton;

  DLOG(INFO) << jsonSkeleton.dump(4);

  file << jsonSkeleton.dump(4);
  file.close();

}

}// namespace ns


#endif