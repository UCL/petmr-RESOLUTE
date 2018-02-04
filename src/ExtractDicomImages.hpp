/*
   ExtractDicomImages.hpp

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

#ifndef _EXTRACTDICOMIMAGES_H_
#define _EXTRACTDICOMIMAGES_H_ 

#include <nlohmann/json.hpp>
#include <boost/filesystem.hpp>
#include <glog/logging.h>

//GDCM includes
#include "gdcmReader.h"

namespace dcm {


const nlohmann::json StudyRecord = R"(
  {
    "StudyUID": " "
  }
  )"_json;

const nlohmann::json SeriesRecord = R"(
  {
    "StudyUID": " ",
    "SeriesUID": " ",
    "SeriesDesc": " ",
    "SeriesNo": 0
  }
)"_json;

const nlohmann::json InstanceRecord = R"(
  {
    "SeriesUID": "",
    "InstanceUID": "",
    "ImageNo": 0,
    "FilePath": ""
  }
)"_json;

typedef std::set<nlohmann::json> InstanceListType, SeriesListType, StudyListType;

bool GetDicomInfo(const boost::filesystem::path srcPath, gdcm::DataSet &ds){

  //Create reader, set filename and check if reading succeeds.
  std::unique_ptr<gdcm::Reader> DICOMreader(new gdcm::Reader);
  DICOMreader->SetFileName(boost::filesystem::canonical(srcPath).string().c_str());

  if (!DICOMreader->Read())
  {
    return false;
  }

  //Extract DICOM header information.
  ds = DICOMreader->GetFile().GetDataSet();

  return true;
}


bool GetTagInfo(const gdcm::DataSet &ds, const gdcm::Tag tag, std::string &dst){

  //Extracts information for a given DICOM tag from a gdcm dataset.
  //Tag contents are returned as a string in dst variable.

  //TODO: Do actual check for valid content.

  //Tries to read the element associated with the tag. If the read fails, the
  //DataElement should have a ByteValue of NULL.

  try {
    std::stringstream inStream;
    inStream.exceptions(std::ios::badbit);
    const gdcm::DataElement& element = ds.GetDataElement(tag);
    if (element.GetByteValue() != NULL) {
      inStream << std::fixed << element.GetValue();
      dst = inStream.str();
    }
    else return false;

  } catch (std::bad_alloc){
    LOG(ERROR) << "DICOMIO::getTagInfo : Cannot read!";
    return false;
  }

  return true;
}


class StudyTree {

public:
  StudyTree(){};
  StudyTree(boost::filesystem::path rootPath){ _rootPath = rootPath; PopulateLists(); };

protected:

  void PopulateLists();
  void AddStudyRecord(const gdcm::DataSet &ds);
  void AddSeriesRecord(const gdcm::DataSet &ds);
  virtual void AddInstanceRecord(const gdcm::DataSet &ds, const boost::filesystem::path pth);

  StudyListType _studyList;
  SeriesListType  _seriesList;
  InstanceListType _instanceList;

  boost::filesystem::path _rootPath;

};

void StudyTree::PopulateLists(){
  _studyList.clear();
  _seriesList.clear();
  _instanceList.clear();

  uint64_t count = 0;

  gdcm::DataSet ds;
  try
  {
    for (auto &entry : boost::make_iterator_range(boost::filesystem::recursive_directory_iterator(_rootPath), {}))
    {
      if (boost::filesystem::is_regular_file(entry.status()))
      {
        //DLOG(INFO) << "Reading " << entry.path();
        //if (inm::IsDICOM(entry.path()))
        //  fileList.push_back(entry.path());
        if (GetDicomInfo(entry.path(), ds)) {
          AddStudyRecord(ds);
          AddSeriesRecord(ds);
          this->AddInstanceRecord(ds, entry.path());
          count++;
        }
      }
    }
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    LOG(ERROR) << "Problem accessing root directory!";
  }
  //End recursion through src directories

  LOG(INFO) << count << " DICOM files found";
  LOG(INFO) << _studyList.size() << " studies found";
  LOG(INFO) << _seriesList.size() << " series found";
  LOG(INFO) << _instanceList.size() << " images found";

  for (auto s : _seriesList){
    DLOG(INFO) << s.dump(4);
  }
  //LOG(INFO) << _studyList.size() << " studies found";

}

void StudyTree::AddStudyRecord(const gdcm::DataSet &ds){

  nlohmann::json study = StudyRecord;

  std::string studyUID;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x00d), studyUID);

  study["StudyUID"] = studyUID;

  _studyList.insert(study);
  //DLOG(INFO) << study.dump(4);

}

void StudyTree::AddSeriesRecord(const gdcm::DataSet &ds){

  nlohmann::json series = SeriesRecord;

  std::string studyUID;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x00d), studyUID);

  std::string seriesUID;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x00e), seriesUID); 

  std::string seriesDesc;
  GetTagInfo(ds,gdcm::Tag(0x0008,0x103e), seriesDesc); 

  std::string seriesNo;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x0011), seriesNo); 

  series["StudyUID"] = studyUID;
  series["SeriesUID"] = seriesUID;
  series["SeriesDesc"] = seriesDesc;
  series["SeriesNo"] = std::stoi(seriesNo);

  _seriesList.insert(series);
  //DLOG(INFO) << study.dump(4);

}

void StudyTree::AddInstanceRecord(const gdcm::DataSet &ds, const boost::filesystem::path pth){

  nlohmann::json instance = InstanceRecord;

  std::string seriesUID;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x00e), seriesUID); 

  std::string instanceUID;
  GetTagInfo(ds,gdcm::Tag(0x0008,0x0018), instanceUID); 

  std::string imageNo;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x0013), imageNo); 

  instance["SeriesUID"] = seriesUID;
  instance["InstanceUID"] = instanceUID;
  instance["ImageNo"] = std::stoi(imageNo);
  instance["FilePath"] = pth.string();

  _instanceList.insert(instance);
  //DLOG(INFO) << instance.dump(4);

}

class UTETree : public StudyTree {
 
public:
    UTETree(boost::filesystem::path rootPath){ _rootPath = rootPath; PopulateLists(); };

protected:
  void AddInstanceRecord(const gdcm::DataSet &ds, const boost::filesystem::path pth) override;  
};

void UTETree::AddInstanceRecord(const gdcm::DataSet &ds, const boost::filesystem::path pth) {

  nlohmann::json instance = InstanceRecord;

  std::string seriesUID;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x00e), seriesUID); 

  std::string instanceUID;
  GetTagInfo(ds,gdcm::Tag(0x0008,0x0018), instanceUID); 

  std::string imageNo;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x0013), imageNo); 

  std::string echoTime;
  GetTagInfo(ds,gdcm::Tag(0x0018,0x0081), echoTime);   

  instance["SeriesUID"] = seriesUID;
  instance["InstanceUID"] = instanceUID;
  instance["ImageNo"] = std::stoi(imageNo);
  instance["FilePath"] = pth.string();
  instance["TE"] = echoTime;

  _instanceList.insert(instance);
  DLOG(INFO) << instance.dump(4);

}


} //namespace dcm


#endif