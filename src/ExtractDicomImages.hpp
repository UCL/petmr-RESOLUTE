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

#include <string>

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

  int GetNoOfStudies(){ return _studyList.size(); };
  std::string GetStudyUID( unsigned int pos );

  int GetNoOfSeries(const std::string studyUID);
  std::vector<std::string> GetSeriesUIDList(const std::string studyUID);

  unsigned int GetNoOfImages(const std::string seriesUID);
  std::vector<nlohmann::json> GetInstanceList(const std::string seriesUID);

  std::vector<boost::filesystem::path> GetSeriesFileList(const std::string seriesUID);


protected:

  void PopulateLists();
  void AddStudyRecord(const gdcm::DataSet &ds);
  void AddSeriesRecord(const gdcm::DataSet &ds);

  nlohmann::json GetSeriesRecord(const std::string seriesUID);

  void GetBasicInstanceInfo(const gdcm::DataSet &ds, nlohmann::json &instance);

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

  LOG_IF(INFO, _studyList.size() == 1) << "\t" << _studyList.size() << " study";
  LOG_IF(INFO, _studyList.size() != 1) << "\t" << _studyList.size() << " studies";

  LOG(INFO) << "\t" << _seriesList.size() << " series";
  LOG(INFO) << "\t" << _instanceList.size() << " images";

  for (auto s : _seriesList){
    DLOG(INFO) << "Series List: " << std::endl << s.dump(4);
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

void StudyTree::GetBasicInstanceInfo(const gdcm::DataSet &ds, nlohmann::json &instance){

  std::string seriesUID;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x00e), seriesUID); 

  std::string instanceUID;
  GetTagInfo(ds,gdcm::Tag(0x0008,0x0018), instanceUID); 

  std::string imageNo;
  GetTagInfo(ds,gdcm::Tag(0x0020,0x0013), imageNo); 

  instance["SeriesUID"] = seriesUID;
  instance["InstanceUID"] = instanceUID;
  instance["ImageNo"] = std::stoi(imageNo);

}

void StudyTree::AddInstanceRecord(const gdcm::DataSet &ds, const boost::filesystem::path pth){

  nlohmann::json instance = InstanceRecord;
  GetBasicInstanceInfo(ds,instance);
  instance["FilePath"] = pth.string();

  _instanceList.insert(instance);
  //DLOG(INFO) << instance.dump(4);

}

int StudyTree::GetNoOfSeries(const std::string studyUID){

  int noFound = 0;

  for (auto const& x : _seriesList){
    if (x["StudyUID"] == studyUID)
      noFound++;
  }

  return noFound;
}

std::string StudyTree::GetStudyUID( unsigned int pos ){

  if ( (pos > _studyList.size()) || (pos == 0)) {
    throw false;
  }

  int target = pos-1;
  int count = 0;

  for (auto const& x : _studyList){
    if (count == target)
      return x["StudyUID"];

    count++;
  }
  return "";

}

std::vector<std::string> StudyTree::GetSeriesUIDList(const std::string studyUID){

  int noOfSeries = GetNoOfSeries(studyUID);

  std::vector<std::string> outList;

  if (noOfSeries == 0){
    LOG(ERROR) << "Study UID: " << studyUID << " not found!";
    throw false; 
  }

  for (auto const& x : _seriesList){
    if (x["StudyUID"] == studyUID)
      outList.push_back(x["SeriesUID"]);
  }

  return outList;
}

unsigned int StudyTree::GetNoOfImages(const std::string seriesUID){

  unsigned int count=0;

  for (auto const& s : _instanceList) {
    if (s["SeriesUID"] == seriesUID)
      count++;
  }
  return count;
}

nlohmann::json StudyTree::GetSeriesRecord(const std::string seriesUID){

  for (auto const& x : _seriesList){
    if (x["SeriesUID"] == seriesUID)
      return x;
  }

  return {};
}

std::vector<nlohmann::json> StudyTree::GetInstanceList(const std::string seriesUID){

  std::vector<nlohmann::json> outList;

  for (auto const& i: _instanceList){
    if (i["SeriesUID"] == seriesUID)
      outList.push_back(i);    
  }

  return outList;
}

std::vector<boost::filesystem::path> StudyTree::GetSeriesFileList(const std::string seriesUID){
  
  std::vector<boost::filesystem::path> outList;

  for (auto const& i: _instanceList){
    if (i["SeriesUID"] == seriesUID){
      std::string path = i["FilePath"];
      //DLOG(INFO) << "\t : " << path; 
      outList.push_back(path);    
    }
  }

  return outList; 
}


class UTETree : public StudyTree {
 
public:
  UTETree(boost::filesystem::path rootPath){ _rootPath = rootPath; PopulateLists(); };

  std::string FindMuMapUID(const std::string studyUID, const std::string tag);
  std::string FindUTEUID(const std::string studyUID, const std::string tag, const std::string TE);

protected:
  void AddInstanceRecord(const gdcm::DataSet &ds, const boost::filesystem::path pth) override; 

  bool CheckSeriesTE(const std::string seriesUID, const std::string TE);

};

void UTETree::AddInstanceRecord(const gdcm::DataSet &ds, const boost::filesystem::path pth) {

  nlohmann::json instance = InstanceRecord;
  GetBasicInstanceInfo(ds,instance);

  std::string echoTime;
  GetTagInfo(ds,gdcm::Tag(0x0018,0x0081), echoTime);   

  instance["FilePath"] = pth.string();
  instance["TE"] = echoTime;

  _instanceList.insert(instance);
  //DLOG(INFO) << instance.dump(4);

}

std::string UTETree::FindMuMapUID(const std::string studyUID, const std::string tag){

  std::vector<std::string> seriesToEval = GetSeriesUIDList(studyUID);

  if (seriesToEval.size() == 0) {
    LOG(ERROR) << "No mu-map found";
    throw false;
    return "";
  }

  for (auto const &s : seriesToEval){
    nlohmann::json seriesRec = GetSeriesRecord(s);
    std::string desc = seriesRec["SeriesDesc"];
    if (desc.find(tag) != std::string::npos){
      LOG(INFO) << "Identified mu-map series: " << seriesRec["SeriesUID"];
      return seriesRec["SeriesUID"];
    }
  }

  LOG(ERROR) << "No mu-map found!";
  throw false;

  return "";
}

bool UTETree::CheckSeriesTE(const std::string seriesUID, const std::string TE){

  std::vector<nlohmann::json> instList = GetInstanceList(seriesUID);

  for (auto const& i : instList){
    if (i["TE"] != TE)
      return false;
  }

  return true;
}

std::string UTETree::FindUTEUID(const std::string studyUID, const std::string tag, const std::string TE){

  std::vector<std::string> seriesToEval = GetSeriesUIDList(studyUID);

  if (seriesToEval.size() == 0) {
    LOG(ERROR) << "No UTE series found";
    throw false;
    return "";
  }

  for (auto const &s : seriesToEval){
    nlohmann::json seriesRec = GetSeriesRecord(s);
    std::string desc = seriesRec["SeriesDesc"];
    if (desc.find(tag) != std::string::npos){

      //IF series has TE
      if (CheckSeriesTE(s,TE)) {
        LOG(INFO) << "Identified UTE series (TE = " << TE << "): " << seriesRec["SeriesUID"];
        return seriesRec["SeriesUID"];
      }
    }
  }

  LOG(ERROR) << "No UTE found!";
  throw false;

  return "";
}


} //namespace dcm


#endif