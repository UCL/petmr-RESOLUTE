# petmr-RESOLUTE 

[![Build Status](https://travis-ci.com/UCL/petmr-RESOLUTE.svg?token=2LGssZ2qj5A4K3LNd3es&branch=master)](https://travis-ci.com/UCL/petmr-RESOLUTE)

Implementation of RESOLUTE pseudo-CT method for the mMR.

## Required packages
- [ANTs](https://github.com/ANTsX/ANTs)
- [ITK](https://itk.org/) (Note. this can be built when compiling ANTs)
- [Boost](http://www.boost.org/) 
- [glog](https://github.com/google/glog)
- [DCMTK](http://dicom.offis.de/). The application `dcmodify` must be available on your path.

## Basic usage
```shell
./resolute -i <DICOMDIR> -j <JSON>
```
where ```<DICOMDIR>``` contains both the UTE and MRAC DICOM data for a given patient and ```<JSON>``` is the JSON configuration file. The application will produce a folder in the output directory specified in the JSON file. The output folder is named using the ```Study UID```, and inside this folder will be a new DICOM series comprising the RESOLUTE MRAC image.
