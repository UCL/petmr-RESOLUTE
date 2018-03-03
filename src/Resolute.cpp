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

#include <itkCastImageFilter.h>
#include <itkImageSeriesWriter.h>
#include <itkNumericSeriesFileNames.h>
#include <itkGDCMSeriesFileNames.h>
#include <gdcmReader.h>
#include <gdcmWriter.h>
#include <gdcmAttribute.h>
#include <gdcmUIDGenerator.h>

#include "EnvironmentInfo.h"
#include "ParamSkeleton.hpp"
#include "ExtractDicomImages.hpp"
#include "Resolute.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using json = nlohmann::json;

typedef itk::Image<float,3> ImageType;

typedef uint16_t OutputPixelType;
typedef itk::Image<OutputPixelType,3> OutputImageType;

const char *const INMUIDROOT = "1.2.826.0.1.3680043.9.6705";

static std::string GetUUID()
{
    //Produces a new INM UID for DICOM data.
    std::string newUID;

    gdcm::UIDGenerator gen;
    gen.SetRoot(INMUIDROOT);
    newUID = (std::string)gen.Generate();

    if (!gen.IsValid(newUID.c_str()))
        return "";

    return newUID;
}

int GetDICOMImageNumber(const fs::path src){

    std::unique_ptr<gdcm::Reader> dicomReader(new gdcm::Reader);
    dicomReader->SetFileName(src.string().c_str());

    if (!dicomReader->Read()) {
      LOG(INFO) << "Unable to read as DICOM file";
      throw false;
    }

    const gdcm::DataSet &ds = dicomReader->GetFile().GetDataSet();

    const gdcm::Tag instanceNumber(0x0020,0x0013);
    const gdcm::DataElement &pdde = ds.GetDataElement(instanceNumber);

    std::stringstream instanceNumberVal;
    instanceNumberVal << (char *) pdde.GetByteValue()->GetPointer();

    DLOG(INFO) << "Image number: " << instanceNumberVal.str();

    int outputNo = -1;
    instanceNumberVal >> outputNo;

    return outputNo;
}

void CreateDICOMSeriesFromMRAC(
    const ImageType::Pointer &img, 
    const std::vector<fs::path> &originalFiles,
    const fs::path destDir){

  typedef itk::CastImageFilter< ImageType, OutputImageType> CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();

  castFilter->SetInput(img);

  try {
    castFilter->Update();
  } catch (itk::ExceptionObject & err) {
      LOG(ERROR) << "Cannot cast from ImageType to OutputImageType!";
      throw false;
  }

  typedef itk::Image<OutputImageType::PixelType, 2> ImageType2D;
  typedef itk::ImageSeriesWriter< OutputImageType, ImageType2D > SeriesWriterType;

  const OutputImageType::RegionType& inputRegion = img->GetLargestPossibleRegion();
  const OutputImageType::SizeType& inputSize = inputRegion.GetSize();
  const OutputImageType::SpacingType& inputSpacing = img->GetSpacing();

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();

  // Generate the file names
  typedef itk::NumericSeriesFileNames OutputNamesGeneratorType;
  OutputNamesGeneratorType::Pointer outputNames = OutputNamesGeneratorType::New();

  fs::path tmpDir = fs::unique_path();

  fs::path dstPath = destDir;
  dstPath /= tmpDir;
  
   //Create our temp destination directory if it doesn't already exist.
  if (!fs::exists(dstPath)){
    try {
      fs::create_directories(dstPath);
    } catch (const fs::filesystem_error &e){
      LOG(ERROR) << "Cannot create tmp folder : " << dstPath;
      throw false;
    }
  }

  fs::path finalPath = destDir;
  finalPath /= "RESOLUTE_MRAC";
  
   //Create out destination directory if it doesn't already exist.
  if (!fs::exists(finalPath)){
    try {
      fs::create_directories(finalPath);
    } catch (const fs::filesystem_error &e){
      LOG(ERROR) << "Cannot create destination folder : " << finalPath;
      throw false;
    }
  }

  std::string seriesFormat = dstPath.string() + "/" + "IM%d.img";
  LOG(INFO) << "SeriesFormat:" << seriesFormat;

  outputNames->SetSeriesFormat (seriesFormat.c_str());
  outputNames->SetStartIndex(1);
  outputNames->SetEndIndex(inputSize[2]);

  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput( castFilter->GetOutput() );
  seriesWriter->SetFileNames( outputNames->GetFileNames() );

  try {
    LOG(INFO) << "Writing raw files to " << dstPath;
    seriesWriter->Update();
  } catch (itk::ExceptionObject & err) {
    LOG(ERROR) << "Cannot write output slices!";
    throw false;
  }

  std::string newSeriesUID = GetUUID();

  std::string exec = "dcmodify";

  for (int x=0; x < originalFiles.size(); x++){
    fs::path outFilePath = dstPath;
    outFilePath /= "mumap-";
    outFilePath += boost::lexical_cast<std::string>(x+1);
    outFilePath += ".dcm";

    try {
      fs::copy(originalFiles[x], outFilePath);
    } catch (const fs::filesystem_error &e){
      LOG(ERROR) << "Cannot copy original MRAC image to " << outFilePath;
      throw false;
    }

    fs::path pctfile = dstPath;
    pctfile /= "IM";
    pctfile += boost::lexical_cast<std::string>(GetDICOMImageNumber(originalFiles[x]));
    pctfile += ".img";

    std::string args = exec + " " + "-nb -if PixelData=\"" + pctfile.string() + "\" \"" + outFilePath.string() + "\"";
    DLOG(INFO) << "args= " << args;
    system( args.c_str() );

    args = exec + " " + "-m 0008,103E=\"RESOLUTE MRAC\"" + " -m 0020,0011=1999 -m 0020,000E=" + newSeriesUID + " -gin \"" + outFilePath.string() + "\"";
    DLOG(INFO) << "args= " << args;
    system( args.c_str() );

    fs::path finalFilePath = finalPath;
    finalFilePath /= newSeriesUID;
    finalFilePath += ".";
    finalFilePath += boost::lexical_cast<std::string>(x+1);
    finalFilePath += ".dcm";

    try {
      fs::rename(outFilePath, finalFilePath);
    } catch (const fs::filesystem_error &e){
      LOG(ERROR) << "Cannot move new MRAC image to " << finalFilePath;
      throw false;
    }    

  }

  //Delete the temp. folder.
  try
  {
    if(fs::exists(dstPath))
      fs::remove_all(dstPath);
  } 
  catch(const fs::filesystem_error & e){
      LOG(WARNING) << "Could not delete temp. directory: " << dstPath;
  }

}

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

  typedef dcm::ReadDicomSeries<ImageType> SeriesReadType;

  typedef ns::ResoluteImageFilter<ImageType,ImageType> ResoluteFilterType;
  ResoluteFilterType::Pointer resoluteFilter = ResoluteFilterType::New();
  resoluteFilter->SetJSONParams(paramFile);
  resoluteFilter->SetOutputDirectory(destRoot);
  resoluteFilter->SetOutputFileExtension(outputType);

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

  const float SIEMENS_VOX_SCALING = 10000.0;
  typedef typename itk::MultiplyImageFilter<ImageType,ImageType> MultiplyFilterType;
  typename MultiplyFilterType::Pointer mult = MultiplyFilterType::New();
  mult->SetInput(resoluteFilter->GetOutput());
  mult->SetConstant(SIEMENS_VOX_SCALING);

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();

  boost::filesystem::path outFileName = destRoot;
  outFileName /= "RESOLUTE-mMR-scaling.nii.gz";
  writer->SetFileName(outFileName.string());
  writer->SetInput(mult->GetOutput());

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not scaled RESOLUTE image!";
    return EXIT_FAILURE;   
  }
  //Modify DICOM data here

  std::vector<fs::path> mracfileNames = tree->GetSeriesFileList(mumapUID);

  fs::path finalDest = destRoot;
  finalDest /= "DICOM";
  CreateDICOMSeriesFromMRAC(mult->GetOutput(), mracfileNames, finalDest);


  //Print total execution time
  std::time_t stopTime = std::time( 0 ) ;
  unsigned int totalTime = stopTime - startTime;
  LOG(INFO) << "Time taken: " << totalTime << " seconds";
  LOG(INFO) << "Ended: " << std::asctime(std::localtime(&stopTime));
  return EXIT_SUCCESS;
}
