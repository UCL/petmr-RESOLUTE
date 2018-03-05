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

## Configuration file

A skeleton JSON file can be created with the command:
```
./resolute --create-json <MYJSONFILE>
```
where ```<MYJSONFILE>``` is an output filename. 
The skeleton will look like this:
```yaml
{
    "MRACSeriesName": "Head_MRAC_PET_UTE_UMAP",
    "UTE1SeriesName": "Head_MRAC_PET_UTE",
    "UTE1TE": "0.07",
    "UTE2SeriesName": "Head_MRAC_PET_UTE",
    "UTE2TE": "2.46",
    "destDir": ".",
    "destExportMethod": "FILE",
    "destFileType": ".nii.gz",
    "logDir": "./logs",
    "regArgs": "3 -m CC[<%%REF%%>,<%%FLOAT%%>,1,4] -i 10x5x2 -o <%%PREFIX%%> -t SyN[0.5] -r Gauss[3,0] -G",
    "regName": "ANTS",
    "regTemplatePath": "",
    "version": "0.0.1"
}
```
Fill in the desired output directory in the variable `destDir`, and the path to the registration template in `regTemplatePath`, e.g.
```yaml
{
    "MRACSeriesName": "Head_MRAC_PET_UTE_UMAP",
    "UTE1SeriesName": "Head_MRAC_PET_UTE",
    "UTE1TE": "0.07",
    "UTE2SeriesName": "Head_MRAC_PET_UTE",
    "UTE2TE": "2.46",
    "destDir": "./OUTPUT",
    "destExportMethod": "FILE",
    "destFileType": ".nii.gz",
    "logDir": "./logs",
    "regName": "ANTS",
    "regTemplatePath": "/path/to/template/mni_icbm152_t1_tal_nlin_sym_09a.nii.gz",
    "regArgs": "3 -m CC[<%%REF%%>,<%%FLOAT%%>,1,4] -i 10x5x2 -o <%%PREFIX%%> -t SyN[0.5] -r Gauss[3,0] -G",
    "version": "0.0.1"
}

```
