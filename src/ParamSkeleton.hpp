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
    double UTE1TR;
    std::string UTE2SeriesName;
    double UTE2TR;

    std::string regName;
    boost::filesystem::path regTemplatePath;
    boost::filesystem::path regMaskPath;
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
        {"UTE1TR", p.UTE1TR},
        {"UTE2SeriesName", p.UTE2SeriesName},
        {"UTE2TR", p.UTE2TR},

        {"regName", p.regName},
        {"regTemplatePath", p.regTemplatePath.string()},
        {"regMaskPath", p.regMaskPath.string()},
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
    p.UTE1TR = j.at("UTE1TR").get<double>();    
    p.UTE2SeriesName = j.at("UTE2SeriesName").get<std::string>();
    p.UTE2TR = j.at("UTE2TR").get<double>();    

    p.regName = j.at("regName").get<std::string>();
    p.regTemplatePath = j.at("regTemplatePath").get<std::string>();
    p.regMaskPath = j.at("regMaskPath").get<std::string>();
    p.regArgs = j.at("regArgs").get<std::string>();

  }

  const params skeleton = {
    VERSION_NO,
    ".nii.gz",
    ".",
    "FILE",

    "./logs",

    "UMAP",
    "UTE",
    0.07,
    "UTE",
    2.46,

    "ANTS",
    "",
    "",
    "",
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

  //Check if registration mask exists
  if (!boost::filesystem::exists(pr.regMaskPath)){
    LOG(ERROR) << "Registration mask " << pr.regMaskPath << " does not exist!";
    throw false;
  }

  //Check if registration mask is a file
  if (!boost::filesystem::is_regular_file(pr.regMaskPath)){
    LOG(ERROR) << "Registration mask " << pr.regMaskPath << " does not appear to be a file!";
    throw false;
  }




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