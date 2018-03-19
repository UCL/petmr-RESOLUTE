/*
   TemplateController.hpp

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

#ifndef _TEMPLATECONTROLLER_HPP_
#define _TEMPLATECONTROLLER_HPP_

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <fstream>
#include <glog/logging.h>
#include <nlohmann/json.hpp>

namespace tc {

enum class ETemplateImages {
  GM, WM, CSF, Brain, Frontal, Mastoid, Nasal, Skull, T1
};

class TemplateController {

public:

  TemplateController(){};
  void SetPath(const boost::filesystem::path &pth);
  boost::filesystem::path GetFilePath(const ETemplateImages e);
  std::string GetFileName(const ETemplateImages e);

protected:

  boost::filesystem::path _rootDir;
  nlohmann::json _jsonManifest;

};

void TemplateController::SetPath(const boost::filesystem::path &pth){

  boost::filesystem::path tempPath;

  if (boost::filesystem::is_regular_file(pth)){
    tempPath = pth;
  }
  else {
    if (boost::filesystem::is_directory(pth)){
      tempPath = pth;
      tempPath /= "manifest.yaml";
    }
  }

  LOG(INFO) << "Reading manifest from: " << tempPath;
  
  try {
    std::ifstream ifs(tempPath.c_str());
    _jsonManifest = nlohmann::json::parse(ifs);
  } catch (...) {
    LOG(ERROR) << "Unable to read manifest!";
    throw false;
  }
  _rootDir = pth.parent_path();

}
boost::filesystem::path TemplateController::GetFilePath(const ETemplateImages e){

  boost::filesystem::path outDir = _rootDir;
  outDir /= GetFileName(e);
  return outDir;

}

std::string TemplateController::GetFileName(const ETemplateImages e){

  switch (e){
    case ETemplateImages::GM: return _jsonManifest["GMReg"].template get<std::string>() ; break;
    case ETemplateImages::WM: return _jsonManifest["WMReg"].template get<std::string>(); break;
    case ETemplateImages::CSF: return _jsonManifest["CSFReg"].template get<std::string>(); break; 
    case ETemplateImages::Brain: return _jsonManifest["brainMask"].template get<std::string>(); break;
    case ETemplateImages::Frontal: return _jsonManifest["frontalReg"].template get<std::string>(); break;
    case ETemplateImages::Mastoid: return _jsonManifest["mastoidReg"].template get<std::string>(); break;
    case ETemplateImages::Nasal: return _jsonManifest["nasalReg"].template get<std::string>(); break;
    case ETemplateImages::Skull: return _jsonManifest["skullReg"].template get<std::string>(); break;
    case ETemplateImages::T1: return _jsonManifest["template"].template get<std::string>(); break;
    default: return "";
  }

  return "";
}

}// end namespace tc

#endif