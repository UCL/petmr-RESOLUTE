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

/*
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageDuplicator.h>
#include <itkExtractImageFilter.h>*/

namespace ns {

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

  void SetMRACImage(const TInputImage* mrac);
  void SetUTEImage1(const TInputImage* ute1);
  void SetUTEImage2(const TInputImage* ute2);
  void SetMaskImage(const TMaskImage* mask);

protected:
  ResoluteImageFilter();
  ~ResoluteImageFilter(){};

  typename TInputImage::ConstPointer GetMRACImage();
  typename TInputImage::ConstPointer GetUTEImage1();
  typename TInputImage::ConstPointer GetUTEImage2();
  typename TMaskImage::ConstPointer GetMaskImage();

  typename HistoImageType::Pointer _histogram;

  struct cluster_coord {
    unsigned int x,y;
  };

  typedef typename std::vector<cluster_coord> CoordListVector;

  typename ResoluteImageFilter::CoordListVector _coords;

  void CalculateHistogram();
  void FindClusterCoords();

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

  writer->SetFileName("histo.mhd");
  writer->SetInput(_histogram);

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write histogram!";
    throw(ex);    
  }

  /*
  typedef itk::ImageDuplicator<OutputHistogramImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( histogramToImageFilter->GetOutput() );

  try {
    duplicator->Update();
    dstImage = duplicator->GetOutput();
  }
  catch (itk::ExceptionObject &ex){
    BOOST_LOG_TRIVIAL(error) << "Could not duplicate histogram image!";
    return false;  
  }*/

}

template< typename TInputImage, typename TMaskImage>
void ResoluteImageFilter<TInputImage, TMaskImage>::FindClusterCoords(){
  
  typedef itk::ThresholdImageFilter< HistoImageType > ThresholdType;
  ThresholdType::Pointer thresh = ThresholdType::New();

  thresh->SetInput( _histogram );
  thresh->ThresholdAbove(5000);
  thresh->SetOutsideValue(0);
  thresh->Update();

  LOG(INFO) << "Starting k-means";

  const unsigned int NUMOFCLASSES=3;
  _coords.clear();

  typedef itk::Statistics::JointDomainImageToListSampleAdaptor< HistoImageType >   AdaptorType;
  typename AdaptorType::Pointer adaptor = AdaptorType::New();
  adaptor->SetImage( thresh->GetOutput() );

  typedef typename AdaptorType::MeasurementVectorType  MeasurementVectorType;

  // Create the K-d tree structure
  typedef itk::Statistics::WeightedCentroidKdTreeGenerator<AdaptorType > TreeGeneratorType;

  typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

  treeGenerator->SetSample( adaptor );
  treeGenerator->SetBucketSize( 8 );
  treeGenerator->Update();

  DLOG(INFO) << "Generated K-d tree";

  typedef typename TreeGeneratorType::KdTreeType TreeType;
  typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;

  typename EstimatorType::Pointer estimator = EstimatorType::New();
  typename EstimatorType::ParametersType initialMeans( 9 );
  
  initialMeans[0] = 0;
  initialMeans[1] = 0;
  initialMeans[2] = 1000;
  initialMeans[3] = 0;
  initialMeans[4] = 0;
  initialMeans[5] = 1000;
  initialMeans[6] = 0;
  initialMeans[7] = 0;
  initialMeans[8] = 0;   

  estimator->SetParameters( initialMeans );
  estimator->SetKdTree( treeGenerator->GetOutput() );
  estimator->SetMaximumIteration( 50000 );
  estimator->SetCentroidPositionChangesThreshold(0);
  estimator->StartOptimization();

  LOG(INFO) << "k-means estimation complete.";

  typename EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
  //LOG(INFO) << "Estimated means: " << estimatedMeans;
  //typename EstimatorType::ParametersType estimatedCentroids;

  LOG(INFO) << "Centre of soft-tissue cluster = (" << (int)estimatedMeans[0]-1 << "," << (int)estimatedMeans[1]-1 << ")";
  float scaleFact1 = 1000.0/estimatedMeans[0];
  float scaleFact2 = 1000.0/estimatedMeans[1];

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

  writer->SetFileName("ute1-test.nii.gz");
  writer->SetInput(mult->GetOutput());

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write histogram!";
    throw(ex);    
  }

  typename TInputImage::ConstPointer ute2 = this->GetUTEImage2();

  mult->SetInput1(ute2);
  mult->SetConstant(scaleFact2);

  writer->SetFileName("ute2-test.nii.gz");
  writer->SetInput(mult->GetOutput());

  try {
    writer->Update();
  } catch (itk::ExceptionObject &ex){
    LOG(ERROR) << "Could not write histogram!";
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

  DLOG(INFO) << "Initialised input images.";
  CalculateHistogram();
  DLOG(INFO) << "Calculated histogram.";

  DLOG(INFO) << "Finding centroid";
  FindClusterCoords();
  DLOG(INFO) << "Centroid found at...";

}

}//namespace ns

#endif