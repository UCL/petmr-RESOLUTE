/*
   Resolute.hpp

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

#ifndef _RESOLUTE_HPP_
#define _RESOLUTE_HPP_

#include <boost/filesystem.hpp>
#include <glog/logging.h>

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include <itkJoinImageFilter.h>
#include <itkImageToHistogramFilter.h>
#include <itkHistogramToIntensityImageFilter.h>
#include <itkMinimumMaximumImageFilter.h>

#include <itkKdTree.h>
#include <itkKdTreeBasedKmeansEstimator.h>
#include <itkScalarImageKmeansImageFilter.h>
#include <itkImageToListSampleAdaptor.h>
#include <itkJointDomainImageToListSampleAdaptor.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkScalarImageKmeansImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryFillholeImageFilter.h>

#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include <itkLogImageFilter.h>

//#include <itkVotingBinaryHoleFillFloodingImageFilter.h>
/*
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageDuplicator.h>
#include <itkExtractImageFilter.h>*/

namespace ns {

inline float GetHUfromR2s(float r){
  //From Ladefoged et al. Figure 1.
  return 1.351e-6 * pow(r,3.0) - 3.617e-3 * pow(r,2.0) + 3.841 * r - 19.46;
}

template< typename TInputImage, typename TMaskImage>
class ResoluteImageFilter : public itk::ImageToImageFilter< TInputImage, TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef ResoluteImageFilter Self;
  typedef itk::ImageToImageFilter< TInputImage, TInputImage > Superclass;
  typedef itk::SmartPointer< Self >  Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ResoluteImageFilter, ImageToImageFilter);

  /** Image related typedefs. */
  typedef TInputImage InputImageType;
  typedef typename TInputImage::ConstPointer  InputImagePointer;
  typedef typename TInputImage::RegionType RegionType;
  typedef typename TInputImage::SizeType   SizeType;
  typedef typename TInputImage::IndexType  IndexType;
  typedef typename TInputImage::PixelType  PixelType;

  /** Mask image related typedefs. */
  typedef TMaskImage  MaskImageType;
  typedef typename TMaskImage::ConstPointer MaskImagePointer;
  typedef typename TMaskImage::RegionType MaskRegionType;
  typedef typename TMaskImage::SizeType   MaskSizeType;
  typedef typename TMaskImage::IndexType  MaskIndexType;
  typedef typename TMaskImage::PixelType  MaskPixelType;

  typedef typename itk::Image<float, 2 > HistoImageType;
  typedef itk::Image< unsigned char, 3 > InternalMaskImageType;


  void SetMRACImage(const TInputImage* mrac);
  void SetUTEImage1(const TInputImage* ute1);
  void SetUTEImage2(const TInputImage* ute2);
  void SetMaskImage(const TMaskImage* mask);

  void SetOutputDirectory(const boost::filesystem::path p){ _dstDir = p;};

protected:
  ResoluteImageFilter();
  ~ResoluteImageFilter(){};

  typename TInputImage::ConstPointer GetMRACImage();
  typename TInputImage::ConstPointer GetUTEImage1();
  typename TInputImage::ConstPointer GetUTEImage2();
  typename TMaskImage::ConstPointer GetMaskImage();

  void MakeAirMask();
  void MakePatientVolumeMask();
  void MakeR2s();

  typename HistoImageType::Pointer _histogram;

  typename TInputImage::Pointer _normUTE1;
  typename TInputImage::Pointer _normUTE2;
  typename TInputImage::Pointer _sumUTE;
  typename TInputImage::Pointer _R2s;


  typename InternalMaskImageType::Pointer _airMask;
  typename InternalMaskImageType::Pointer _patVolMask;

  boost::filesystem::path _dstDir;

  struct cluster_coord {
    unsigned int x,y;
  };

  typedef typename std::vector<cluster_coord> CoordListVector;

  //typename ResoluteImageFilter::CoordListVector _coords;
  cluster_coord _coords;

  void CalculateHistogram();
  void GetKMeansMask(const HistoImageType::Pointer &h, HistoImageType::Pointer &outputImage);
  void FindClusterCoords();
  void NormaliseUTE();

  virtual void GenerateData() override;

  ResoluteImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

};

template< typename TInputImage, typename TMaskImage>
ResoluteImageFilter<TInputImage,TMaskImage>
::ResoluteImageFilter()
{
  this->SetNumberOfRequiredInputs(4);
}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::SetMRACImage(const TInputImage* image)
{
  this->SetNthInput(0, const_cast<TInputImage*>(image));
}
 
template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::SetUTEImage1(const TInputImage* image)
{
  this->SetNthInput(1, const_cast<TInputImage*>(image));
}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::SetUTEImage2(const TInputImage* image)
{
  this->SetNthInput(2, const_cast<TInputImage*>(image));
}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::SetMaskImage(const TMaskImage* image)
{
  this->SetNthInput(3, const_cast<TMaskImage*>(image));
}

template< typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer ResoluteImageFilter<TInputImage, TMaskImage>::GetMRACImage()
{
  return static_cast< const TInputImage * >
         ( this->ProcessObject::GetInput(0) );
}

template< typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer ResoluteImageFilter<TInputImage, TMaskImage>::GetUTEImage1()
{
  return static_cast< const TInputImage * >
         ( this->ProcessObject::GetInput(1) );
}

template< typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer ResoluteImageFilter<TInputImage, TMaskImage>::GetUTEImage2()
{
  return static_cast< const TInputImage * >
         ( this->ProcessObject::GetInput(2) );
}

template< typename TInputImage, typename TMaskImage>
typename TMaskImage::ConstPointer ResoluteImageFilter<TInputImage, TMaskImage>::GetMaskImage()
{
  return static_cast< const TMaskImage * >
         ( this->ProcessObject::GetInput(3) );
}


template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::CalculateHistogram()
{

  typedef itk::JoinImageFilter<TInputImage, TInputImage> JoinFilterType;
  typedef typename JoinFilterType::OutputImageType VectorImageType;
  typedef itk::Statistics::ImageToHistogramFilter< VectorImageType > HistogramFilterType;
  typedef itk::MinimumMaximumImageFilter<TInputImage> MinMaxFilterType;

  typename TInputImage::ConstPointer ute1 = this->GetUTEImage1();
  typename TInputImage::ConstPointer ute2 = this->GetUTEImage2();

  typename JoinFilterType::Pointer joinFilter = JoinFilterType::New();

  joinFilter->SetInput1( ute1 );
  joinFilter->SetInput2( ute2 );

  try {
    joinFilter->Update();
  }
  catch (itk::ExceptionObject &e) {
    LOG(ERROR) << e;
    throw(e);
  }

  typename HistogramFilterType::Pointer histoFilter = HistogramFilterType::New();
  typename MinMaxFilterType::Pointer minmaxFilter = MinMaxFilterType::New();

  histoFilter->SetInput( joinFilter->GetOutput() );
  histoFilter->SetMarginalScale(1);

  typedef typename HistogramFilterType::HistogramSizeType HistogramSizeType;
  HistogramSizeType size(2);

  typedef typename HistogramFilterType::HistogramMeasurementVectorType HistogramMeasurementVectorType;
  HistogramMeasurementVectorType binMin(2);
  HistogramMeasurementVectorType binMax(2);

  minmaxFilter->SetInput( ute1 );
  try {
    minmaxFilter->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not calculate min/max for UTE 1!";
    throw(ex);    
  }


  DLOG(INFO) << "UTE 1: min = " << minmaxFilter->GetMinimum();

  binMin[0]=minmaxFilter->GetMinimum();
  binMax[0]=minmaxFilter->GetMaximum();

  minmaxFilter->SetInput( ute2 );
  try {
    minmaxFilter->Update(); 
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not calculate min/max for UTE 2!";
    throw(ex);    
  }

  binMin[1]=minmaxFilter->GetMinimum();
  binMax[1]=minmaxFilter->GetMaximum();

  DLOG(INFO) << "UTE 2: min = " << minmaxFilter->GetMinimum();

  size[0] = static_cast<unsigned int>( binMax[0] - binMin[0] + 1 );
  size[1] = static_cast<unsigned int>( binMax[1] - binMin[1] + 1 );

  histoFilter->SetHistogramSize( size );
  histoFilter->SetHistogramBinMinimum( binMin );
  histoFilter->SetHistogramBinMaximum( binMax );

  try {
    histoFilter->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not calculate joint histogram!";
    throw(ex);    
  }

  typedef typename HistogramFilterType::HistogramType HistogramType;
  const HistogramType *histogram = histoFilter->GetOutput();

  typedef typename itk::HistogramToIntensityImageFilter< HistogramType, HistoImageType > HistogramToImageFilterType;
  typename HistogramToImageFilterType::Pointer histogramToImageFilter = HistogramToImageFilterType::New();

  histogramToImageFilter->SetInput( histogram );

  try {
    histogramToImageFilter->Update();
  }
  catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not convert histogram to image!";
    throw(ex);    
  }

  _histogram = histogramToImageFilter->GetOutput();

  typedef itk::ImageFileWriter<HistoImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  boost::filesystem::path outFileName = _dstDir;
  outFileName /= "histogram.mhd";
  writer->SetFileName(outFileName.string());
  writer->SetInput(_histogram);

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write histogram!";
    throw(ex);    
  }

}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::MakeAirMask(){

  //Takes summed UTE images and creates binary mask of all voxels < 600 (air).
  typedef itk::BinaryThresholdImageFilter <TInputImage, InternalMaskImageType> BinThresholdImageFilterType;
  typename BinThresholdImageFilterType::Pointer binFilter = BinThresholdImageFilterType::New();

  binFilter->SetInput( _sumUTE );
  binFilter->SetLowerThreshold( 0 );
  binFilter->SetUpperThreshold( 600 );
  binFilter->SetInsideValue( 1 );
  binFilter->Update();

  typedef itk::ImageDuplicator<InternalMaskImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(binFilter->GetOutput());
  duplicator->Update();
  _airMask = duplicator->GetOutput();

  typedef itk::ImageFileWriter<InternalMaskImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();

  boost::filesystem::path outFileName = _dstDir;
  outFileName /= "air.nii.gz";
  writer->SetFileName(outFileName.string());
  writer->SetInput(_airMask);

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write air mask!";
    throw(ex);    
  }

}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::MakePatientVolumeMask(){

  //Creates the head mask from MRAC and UTE.

  //MRAC after region growing = 1
  //snUTE > 1000 = 1
  //Add both masks and binarize.
  const int RADIUS = 11;

  typename TInputImage::ConstPointer mrac = this->GetMRACImage();

  typedef itk::BinaryThresholdImageFilter <TInputImage, InternalMaskImageType> BinThresholdImageFilterType;
  typename BinThresholdImageFilterType::Pointer binFilter = BinThresholdImageFilterType::New();
  typename BinThresholdImageFilterType::Pointer binFilter2 = BinThresholdImageFilterType::New();

  binFilter->SetInput( mrac );
  binFilter->SetLowerThreshold( 1 );
  binFilter->SetInsideValue( 1 );

  typedef itk::BinaryBallStructuringElement<InternalMaskImageType::PixelType, InternalMaskImageType::ImageDimension>
              StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(RADIUS);
  structuringElement.CreateStructuringElement();
  typedef itk::BinaryMorphologicalClosingImageFilter <InternalMaskImageType, InternalMaskImageType, StructuringElementType>
          BinaryMorphologicalClosingImageFilterType;
  BinaryMorphologicalClosingImageFilterType::Pointer closingFilter
          = BinaryMorphologicalClosingImageFilterType::New();
  closingFilter->SetInput(binFilter->GetOutput());
  closingFilter->SetForegroundValue(1);
  closingFilter->SetKernel(structuringElement);
  closingFilter->Update();

  /*
  typedef itk::BinaryFillholeImageFilter< CompImageType > FillHoleFilterType;
  typename FillHoleFilterType::Pointer fillFilter = FillHoleFilterType::New();
  fillFilter->SetInput(binFilter->GetOutput());
  fillFilter->SetForegroundValue( 0 );
  fillFilter->UpdateLargestPossibleRegion();
  fillFilter->Update();
  */

  //Takes summed UTE images and creates binary mask of all voxels > 1000 (soft-tissue).
  binFilter2->SetInput( _sumUTE );
  binFilter2->SetLowerThreshold( 1000 );
  binFilter2->SetInsideValue( 1 );

  typedef itk::AddImageFilter<InternalMaskImageType,InternalMaskImageType> AddFilterType;
  typename AddFilterType::Pointer addFilter = AddFilterType::New();

  addFilter->SetInput1(closingFilter->GetOutput());
  addFilter->SetInput2(binFilter2->GetOutput());
  addFilter->Update();

  typedef itk::ConnectedComponentImageFilter <InternalMaskImageType, InternalMaskImageType >
    ConnectedComponentImageFilterType;
 
  ConnectedComponentImageFilterType::Pointer connected =
    ConnectedComponentImageFilterType::New ();
  connected->SetInput(addFilter->GetOutput());
  connected->Update();
 
  LOG(INFO) << "Number of objects: " << connected->GetObjectCount() << std::endl;

  typedef itk::ImageFileWriter<InternalMaskImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();

  boost::filesystem::path outFileName = _dstDir;
  outFileName /= "patient_vol.nii.gz";
  writer->SetFileName(outFileName.string());
  writer->SetInput(connected->GetOutput());

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write patient volume mask!";
    throw(ex);    
  }
  
  typedef itk::ImageDuplicator<InternalMaskImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(addFilter->GetOutput());
  duplicator->Update();
  _patVolMask = duplicator->GetOutput();
  


}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::MakeR2s(){

  float dUTE = (2.46 - 0.07)/1000.0;

  typename TInputImage::ConstPointer ute1 = this->GetUTEImage1();
  typename TInputImage::ConstPointer ute2 = this->GetUTEImage2();

  typedef itk::LogImageFilter<TInputImage,TInputImage> LogFilterType;
  //typedef itk::ImageDuplicator<InternalMaskImageType> DuplicatorType;

  typename LogFilterType::Pointer l1 = LogFilterType::New();
  typename LogFilterType::Pointer l2 = LogFilterType::New();

  l1->SetInput(ute1);
  l2->SetInput(ute2);

  typedef itk::SubtractImageFilter<TInputImage,TInputImage> SubtractFilterType;
  typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();

  subFilter->SetInput1(l1->GetOutput());
  subFilter->SetInput2(l2->GetOutput());

  typedef itk::DivideImageFilter<TInputImage,TInputImage, TInputImage> DivideFilterType;
  typename DivideFilterType::Pointer divideFilter = DivideFilterType::New();

  divideFilter->SetInput(subFilter->GetOutput());
  divideFilter->SetConstant(dUTE);

  typedef itk::MaskImageFilter< TInputImage, InternalMaskImageType, TInputImage > MaskFilterType;
  typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  maskFilter->SetInput(divideFilter->GetOutput());
  maskFilter->SetOutsideValue( 0 ); 
  maskFilter->SetMaskImage( _patVolMask );

  typedef itk::ImageFileWriter<TInputImage> WriterType;
  typename WriterType::Pointer writer = WriterType::New();

  boost::filesystem::path outFileName = _dstDir;
  outFileName /= "R2s.nii.gz";
  writer->SetFileName(outFileName.string());
  writer->SetInput(maskFilter->GetOutput());

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write R2* image!";
    throw(ex);    
  }
  
  typedef itk::ImageDuplicator<TInputImage> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(divideFilter->GetOutput());
  duplicator->Update();
  _R2s = duplicator->GetOutput();



}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::GetKMeansMask(const HistoImageType::Pointer &h, HistoImageType::Pointer &outputImage ){

  LOG(INFO) << "Starting k-means";

  const unsigned int NUMOFCLASSES=2;
  _coords.x = 0;
  _coords.y = 0;

  typedef itk::Statistics::JointDomainImageToListSampleAdaptor< HistoImageType >   AdaptorType;
  typename AdaptorType::Pointer adaptor = AdaptorType::New();
  adaptor->SetImage( h );

  typedef typename AdaptorType::MeasurementVectorType  MeasurementVectorType;

  // Create the K-d tree structure
  typedef itk::Statistics::WeightedCentroidKdTreeGenerator< AdaptorType > TreeGeneratorType;

  typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

  treeGenerator->SetSample( adaptor );
  treeGenerator->SetBucketSize( 16 );
  treeGenerator->Update();

  DLOG(INFO) << "Generated K-d tree";

  typedef typename TreeGeneratorType::KdTreeType TreeType;
  typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;

  typename EstimatorType::Pointer estimator = EstimatorType::New();
  typename EstimatorType::ParametersType initialMeans( NUMOFCLASSES*3 );
  
  initialMeans[0] = 100;
  initialMeans[1] = 100;
  initialMeans[2] = 1000;
  initialMeans[3] = 100;
  initialMeans[4] = 100;
  initialMeans[5] = 0;
  /*initialMeans[6] = 100;
  initialMeans[7] = 100;
  initialMeans[8] = 0; */

  estimator->SetParameters( initialMeans );
  estimator->SetKdTree( treeGenerator->GetOutput() );
  estimator->SetMaximumIteration( 500 );
  estimator->SetCentroidPositionChangesThreshold(0);
  estimator->StartOptimization();

  LOG(INFO) << "k-means estimation complete.";

  typename EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
  //LOG(INFO) << "Estimated means: " << estimatedMeans;
  //typename EstimatorType::ParametersType estimatedCentroids;

  LOG(INFO) << "Centre of soft-tissue cluster = (" << (int)estimatedMeans[0]-1 << "," << (int)estimatedMeans[1]-1 << ")";
  _coords.x = (int)estimatedMeans[0]-1;
  _coords.y = (int)estimatedMeans[1]-1;

  //Test histo start
  typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
  typedef MembershipFunctionType::Pointer MembershipFunctionPointer;
  typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
  typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

  typedef itk::Statistics::SampleClassifierFilter< AdaptorType > ClassifierType;
  ClassifierType::Pointer classifier = ClassifierType::New();
 
  classifier->SetDecisionRule(decisionRule);
  classifier->SetInput( adaptor );
  classifier->SetNumberOfClasses( NUMOFCLASSES );
 
  typedef ClassifierType::ClassLabelVectorObjectType               ClassLabelVectorObjectType;
  typedef ClassifierType::ClassLabelVectorType                     ClassLabelVectorType;
  typedef ClassifierType::MembershipFunctionVectorObjectType       MembershipFunctionVectorObjectType;
  typedef ClassifierType::MembershipFunctionVectorType             MembershipFunctionVectorType;
 
  ClassLabelVectorObjectType::Pointer  classLabelsObject = ClassLabelVectorObjectType::New();
  classifier->SetClassLabels( classLabelsObject );
 
  ClassLabelVectorType &  classLabelsVector = classLabelsObject->Get();
  classLabelsVector.push_back( 1 );
  classLabelsVector.push_back( 0 );
  //classLabelsVector.push_back( 300 );

  MembershipFunctionVectorObjectType::Pointer membershipFunctionsObject =
    MembershipFunctionVectorObjectType::New();
  classifier->SetMembershipFunctions( membershipFunctionsObject );
 
  MembershipFunctionVectorType &  membershipFunctionsVector = membershipFunctionsObject->Get();
 
  MembershipFunctionType::CentroidType origin( adaptor->GetMeasurementVectorSize() );
  int index = 0;
  for ( unsigned int i = 0 ; i < NUMOFCLASSES ; i++ )
    {
    MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();
    for ( unsigned int j = 0 ; j < adaptor->GetMeasurementVectorSize(); j++ )
      {
      origin[j] = estimatedMeans[index++];
      }
    membershipFunction->SetCentroid( origin );
    membershipFunctionsVector.push_back( membershipFunction.GetPointer() );
    }
 
  classifier->Update();

  // Prepare the output image
  //HistoImageType::Pointer outputImage = HistoImageType::New();
  //outputImage->SetRegions(h->GetLargestPossibleRegion());
  //outputImage->Allocate();
  typedef itk::ImageDuplicator<HistoImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(h);
  duplicator->Update();
  outputImage = duplicator->GetOutput();
  outputImage->FillBuffer(0);
 
  const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
  ClassifierType::MembershipSampleType::ConstIterator iter = membershipSample->Begin();

  // Setup the output image iterator - this is automatically synchronized with the membership iterator since the sample is an adaptor
  itk::ImageRegionIteratorWithIndex<HistoImageType> outputIterator(outputImage,outputImage->GetLargestPossibleRegion());
  outputIterator.GoToBegin();

  while ( iter != membershipSample->End() )
  {
    //DLOG(INFO) << "measurement vector = " << iter.GetMeasurementVector()
    //          << "class label = " << iter.GetClassLabel();

    int classLabel = iter.GetClassLabel();
    outputIterator.Set(classLabel);
    ++iter;
    ++outputIterator;
  }

  typedef itk::ImageFileWriter< HistoImageType  > HistoWriterType;
  HistoWriterType::Pointer outputWriter = HistoWriterType::New();

  boost::filesystem::path outFileName = _dstDir;
  outFileName /= "k-means.mhd";
  outputWriter->SetFileName(outFileName.string());
  outputWriter->SetInput(outputImage);
  outputWriter->Update();
  //Test histo end

}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::FindClusterCoords(){
  
  typedef itk::ThresholdImageFilter< HistoImageType > ThresholdType;
  ThresholdType::Pointer thresh = ThresholdType::New();

  thresh->SetInput( _histogram );
  thresh->ThresholdAbove(2000);
  thresh->SetOutsideValue(2000);
  thresh->Update();

  HistoImageType::Pointer histoMaskImage = HistoImageType::New();

  GetKMeansMask(thresh->GetOutput(), histoMaskImage);

  typedef itk::Image< unsigned short, 2 > CompImageType;
  typedef itk::ConnectedComponentImageFilter <HistoImageType, CompImageType >
    ConnectedComponentImageFilterType;
 
  ConnectedComponentImageFilterType::Pointer connected =
    ConnectedComponentImageFilterType::New ();
  connected->SetInput(histoMaskImage);
  connected->Update();
 
  LOG(INFO) << "Number of objects: " << connected->GetObjectCount() << std::endl;

  typedef itk::RelabelComponentImageFilter<CompImageType, CompImageType> RelabelFilterType;
  RelabelFilterType::Pointer relabelFilter =
    RelabelFilterType::New();
  relabelFilter->SetInput(connected->GetOutput());

  relabelFilter->Update();

  typedef itk::BinaryThresholdImageFilter <CompImageType, CompImageType> BinThresholdImageFilterType;
  typename BinThresholdImageFilterType::Pointer binFilter = BinThresholdImageFilterType::New();

  binFilter->SetInput( relabelFilter->GetOutput() );
  binFilter->SetLowerThreshold( 1 );
  binFilter->SetUpperThreshold( 1 );
  binFilter->SetInsideValue( 1 );
  binFilter->Update();

  typedef itk::MaskImageFilter< HistoImageType, CompImageType, HistoImageType > MaskFilterType;
  typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  maskFilter->SetInput(_histogram);
  maskFilter->SetOutsideValue( 0 ); 
  maskFilter->SetMaskImage( binFilter->GetOutput() );
  maskFilter->Update();

  GetKMeansMask(maskFilter->GetOutput(), histoMaskImage);

}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::NormaliseUTE(){

  float scaleFact1 = 1000.0/_coords.x;
  float scaleFact2 = 1000.0/_coords.y;

  if ( (scaleFact1 > 1.0) && (scaleFact2 > 1.0) ){
    LOG(INFO) <<  "\tUTE1 * " <<  scaleFact1;
    LOG(INFO) <<  "\tUTE2 * " <<  scaleFact2;
  }

  typename TInputImage::ConstPointer ute1 = this->GetUTEImage1();

  typedef typename itk::MultiplyImageFilter<TInputImage,TInputImage> MultiplyType;
  typename MultiplyType::Pointer mult = MultiplyType::New();
  mult->SetInput1(ute1);
  mult->SetConstant(scaleFact1);

  typedef itk::ImageFileWriter<TInputImage> WriterType;
  typename WriterType::Pointer writer = WriterType::New();

  boost::filesystem::path outFileName = _dstDir;
  outFileName /= "ute1.nii.gz";
  writer->SetFileName(outFileName.string());
  writer->SetInput(mult->GetOutput());

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write UTE1!";
    throw(ex);    
  }

  typedef itk::ImageDuplicator<TInputImage> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(mult->GetOutput());
  duplicator->Update();

  _normUTE1 = duplicator->GetOutput();

  typename TInputImage::ConstPointer ute2 = this->GetUTEImage2();

  mult->SetInput1(ute2);
  mult->SetConstant(scaleFact2);

  outFileName = _dstDir;
  outFileName /= "ute2.nii.gz";
  writer->SetFileName(outFileName.string());
  writer->SetInput(mult->GetOutput());

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write UTE2!";
    throw(ex);    
  }

  duplicator->SetInputImage(mult->GetOutput());
  duplicator->Update();

  _normUTE2 = duplicator->GetOutput();

  typedef itk::AddImageFilter<TInputImage,TInputImage> AddFilterType;
  typename AddFilterType::Pointer addFilter = AddFilterType::New();

  addFilter->SetInput1(_normUTE1);
  addFilter->SetInput2(_normUTE2);
  addFilter->Update();

  duplicator->SetInputImage(addFilter->GetOutput());
  duplicator->Update();

  _sumUTE = duplicator->GetOutput();

  outFileName = _dstDir;
  outFileName /= "snUTE-test.nii.gz";
  writer->SetFileName(outFileName.string());
  writer->SetInput(_sumUTE);

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write snUTE!";
    throw(ex);    
  }  

}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::GenerateData()
{
  typename TInputImage::ConstPointer mrac = this->GetMRACImage();
  typename TInputImage::ConstPointer ute1 = this->GetUTEImage1();
  typename TInputImage::ConstPointer ute2 = this->GetUTEImage2();
  typename TMaskImage::ConstPointer mask = this->GetMaskImage();
 
  typename TInputImage::Pointer output = this->GetOutput();
  output->SetRegions(mrac->GetLargestPossibleRegion());
  output->Allocate();

  LOG(INFO) << "Initialised input images.";
  CalculateHistogram();
  LOG(INFO) << "Calculated histogram.";

  LOG(INFO) << "Finding centroid";
  FindClusterCoords();
  LOG(INFO) << "Centroid found at...";

  LOG(INFO) << "Normalising UTE";
  NormaliseUTE();
  LOG(INFO) << "Normalisation complete.";

  LOG(INFO) << "Calculating air mask";
  MakeAirMask();
  LOG(INFO) << "Air mask calculation complete.";

  LOG(INFO) << "Calculating patient volume";
  MakePatientVolumeMask();
  LOG(INFO) << "Patient volume calculation complete.";

  LOG(INFO) << "Calculating R2*";
  MakeR2s();
  LOG(INFO) << "R2* calculation complete.";

}

}//namespace ns

#endif