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

   This program performs RESOLUTE for mMR data.
 */

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/program_options.hpp>
#include <glog/logging.h>
#include <nlohmann/json.hpp>

#include "EnvironmentInfo.h"
#include "ParamSkeleton.hpp"
#include "ExtractDicomImages.hpp"
#include "ExtractMaskImages.hpp"
#include "Resolute.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using json = nlohmann::json;

int main(int argc, char **argv)
{

  const char* APP_NAME = "resolute";

  std::string inputDirectoryPath;
  std::string logPath;
  std::string jsonFile;
  std::string outputDirectory;
  std::string prefixName;

  //Set-up command line options
  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print help information")
    ("version","Print version number")
    //("verbose,v", "Be verbose")
    ("input,i", po::value<std::string>(&inputDirectoryPath), "Input DICOMDIR")
    ("log,l", po::value<std::string>(&logPath), "Write log file")
    ("json,j", po::value<std::string>(&jsonFile),  "Use JSON config file")
    ("create-json", po::value<std::string>(&jsonFile),  "Write config JSON skeleton");


  //Evaluate command line options
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc),
      vm); // can throw

    /** --help option
    */
    if (vm.count("help")) {
      std::cout << APP_NAME << std::endl
        << desc << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("version") ) {
      std::cout << APP_NAME << " : v" << VERSION_NO << std::endl;
      return EXIT_SUCCESS;
    }

    if ( (!vm.count("input")) && (!vm.count("create-json")) ){
      std::cout << APP_NAME << std::endl
        << desc << std::endl;
      return EXIT_SUCCESS;   
    }

    po::notify(vm); // throws on error

  } catch (po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << desc << std::endl;
    return EXIT_FAILURE;
  }

  if (vm.count("create-json") ) {
    try {
      ns::WriteJSONSkeleton(jsonFile);
    } catch (bool) {
      std::cerr << "ERROR: Aborting!" << std::endl;
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  std::ifstream ifs(jsonFile);
  json paramFile = json::parse(ifs);

  //Pretty coloured logging (if supported)
  FLAGS_colorlogtostderr = 1;
  FLAGS_alsologtostderr = 1;

  if (vm.count("log")){
    paramFile["logDir"] = logPath;
  }

  DLOG(INFO) << paramFile;
  
  try {
    ns::ValidateJSON(paramFile);
  } catch(bool) {
    LOG(ERROR) << "Invalid JSON file!";
    LOG(ERROR) << "Aborting!";
    return EXIT_FAILURE;
  }

  //Configure logging
  fs::path newLogPath = fs::complete(paramFile["logDir"].get<std::string>());
  newLogPath /= APP_NAME;
  newLogPath += "-";

  google::InitGoogleLogging(argv[0]);
  google::SetLogDestination(google::INFO, newLogPath.string().c_str());

  std::time_t startTime = std::time( 0 ) ;

  //Application starts here
  LOG(INFO) << "Started: " << std::asctime(std::localtime(&startTime));
  LOG(INFO) << "Running '" << APP_NAME << "' version: " << VERSION_NO;
  LOG(INFO) << "Log path = " << newLogPath;
  LOG(INFO) << "Read JSON parameter file: " << jsonFile << std::endl << paramFile.dump(4);

  fs::path srcPath = inputDirectoryPath;
  
  //Check if input path exists
  if (! fs::exists( srcPath ) )
  {
    LOG(ERROR) << "Input path: " << srcPath << " does not exist!";
    return EXIT_FAILURE;
  }

  //Check if it is a directory.
  if (! fs::is_directory( srcPath ) )
  {
    LOG(ERROR) << srcPath << " does not appear to be a directory!";
    return EXIT_FAILURE;
  }

  LOG(INFO) << "Input directory: " << fs::complete(srcPath);

  //Create DICOM UTE search object.
  std::unique_ptr<dcm::UTETree> tree(new dcm::UTETree(srcPath));

  //Total number of series found for first UID.
  std::string studyUID = tree->GetStudyUID(1);
  LOG(INFO) << "No. series in tree: " << tree->GetNoOfSeries(studyUID);

  //Get all Series UIDs associated with study
  //std::vector<std::string> foundSeriesUIDs = tree->GetSeriesUIDList(tree->GetStudyUID(1));

  std::string mumapUID;
  std::string ute1UID;
  std::string ute2UID;

  try {
    //Find the mu-map
    mumapUID = tree->FindMuMapUID(studyUID, paramFile["MRACSeriesName"]);
  }
  catch (bool){
    LOG(ERROR) << "Could not find mu-map series with description \'" << paramFile["MRACSeriesName"] << "\'";
    LOG(ERROR) << "Aborting!";
    return EXIT_FAILURE;
  }

  try {
    //Find UTE 1
    ute1UID = tree->FindUTEUID(studyUID, paramFile["UTE1SeriesName"], paramFile["UTE1TE"].get<std::string>());
  }
  catch (bool){
    LOG(ERROR) << "Could not find UTE1 series with description \'" << paramFile["UTE1SeriesName"] << "\' and TE = " << paramFile["UTE1TE"];
    LOG(ERROR) << "Aborting!";
    return EXIT_FAILURE;
  }  

  try {
    //Find UTE 2
    ute2UID = tree->FindUTEUID(studyUID, paramFile["UTE2SeriesName"], paramFile["UTE2TE"].get<std::string>());
  }
  catch (bool){
    LOG(ERROR) << "Could not find UTE2 series with description \'" << paramFile["UTE2SeriesName"] << "\' and TE = " << paramFile["UTE2TE"];
    LOG(ERROR) << "Aborting!";
    return EXIT_FAILURE;
  }  

  fs::path destRoot = paramFile["destDir"].get<std::string>();
  destRoot /= studyUID;

  //Create out destination directory if it doesn't already exist.
  if (!fs::exists(destRoot)){
    try {
      fs::create_directories(destRoot);
    } catch (const fs::filesystem_error &e){
      LOG(ERROR) << " cannot create destination folder : " << destRoot;
      return EXIT_FAILURE;
    }
  }

  std::string outputType = paramFile["destFileType"];

  std::vector<std::string> srcSeriesUIDs = {mumapUID, ute1UID, ute2UID};

  typedef itk::Image<float,3> ImageType;
  typedef dcm::ReadDicomSeries<ImageType> SeriesReadType;

  //Read mu-map and both UTEs and write file.

  /*
  for (auto const& imgUID: srcSeriesUIDs) {
    std::vector<fs::path> fNames = tree->GetSeriesFileList(imgUID);
    std::unique_ptr<SeriesReadType> dcm(new SeriesReadType(fNames));

    try {
      dcm->Read();
    } catch(bool){
      LOG(ERROR) << "Could not read series: " << imgUID;
      LOG(ERROR) << "Aborting!";
      return EXIT_FAILURE;      
    }

    fs::path outFilePath = fs::canonical(destRoot);
    outFilePath /= imgUID;
    outFilePath += outputType;

    try {
      dcm->Write(outFilePath);
    } catch(bool){
      LOG(ERROR) << "Could not write series: " << imgUID << " to " << outFilePath;
      LOG(ERROR) << "Aborting!";
      return EXIT_FAILURE;      
    }

  }*/

  typedef ns::ResoluteImageFilter<ImageType,ImageType> ResoluteFilterType;
  ResoluteFilterType::Pointer resoluteFilter = ResoluteFilterType::New();
  resoluteFilter->SetJSONParams(paramFile);
  resoluteFilter->SetOutputDirectory(destRoot);

  std::vector<fs::path> fNames = tree->GetSeriesFileList(mumapUID);
  std::unique_ptr<SeriesReadType> dcm(new SeriesReadType(fNames));

  try {
    dcm->Read();
  } catch(bool){
      LOG(ERROR) << "Could not read mu-map series: " << mumapUID;
      LOG(ERROR) << "Aborting!";
      return EXIT_FAILURE;      
  }

  resoluteFilter->SetMRACImage(dcm->GetOutput());

  fNames = tree->GetSeriesFileList(ute1UID);
  dcm.reset(new SeriesReadType(fNames));

  try {
    dcm->Read();
  } catch(bool){
      LOG(ERROR) << "Could not read UTE1 series: " << mumapUID;
      LOG(ERROR) << "Aborting!";
      return EXIT_FAILURE;      
  }

  resoluteFilter->SetUTEImage1(dcm->GetOutput());

  fNames = tree->GetSeriesFileList(ute2UID);
  dcm.reset(new SeriesReadType(fNames));

  try {
    dcm->Read();
  } catch(bool){
      LOG(ERROR) << "Could not read UTE2 series: " << mumapUID;
      LOG(ERROR) << "Aborting!";
      return EXIT_FAILURE;      
  }

  resoluteFilter->SetUTEImage2(dcm->GetOutput());
  resoluteFilter->SetMaskImage(dcm->GetOutput());


  try {
    resoluteFilter->Update();
  } catch (itk::ExceptionObject &e) {
    LOG(ERROR) << e;
    LOG(ERROR) << "Failed to apply RESOLUTE filter!";
    LOG(ERROR) << "Aborting!";
    return EXIT_FAILURE;    
  }

  //Modify DICOM data here


  //Print total execution time
  std::time_t stopTime = std::time( 0 ) ;
  unsigned int totalTime = stopTime - startTime;
  LOG(INFO) << "Time taken: " << totalTime << " seconds";
  LOG(INFO) << "Ended: " << std::asctime(std::localtime(&stopTime));
  return EXIT_SUCCESS;
}
