

#fields: Tools for Spatial Data

This repository contains some supplemental material in this top level and see the subdirectory **fields** for the head of the standard R package. 
The most current package on CRAN  is listed here as *fields_VERSION.tar.gz* . At the time of writing the version is *9.7* . To get precompiled binaries of the stable versions, however,  you can use  [CRAN](https://cran.r-project.org/web/packages/fields).

To create a possibly new version from this repository download the 
fields *subdirectory* and in UNIX. 

```
 R CMD build --force fields
```
#About fields

This is an R package 
 for curve, surface and function fitting with an emphasis
 on splines, spatial data and spatial statistics. The major methods
 include cubic, and thin plate splines, Kriging, and compactly supported
 covariance functions for large data sets. The splines and Kriging methods are
 supported by functions that can determine the smoothing parameter
 (nugget and sill variance) and other covariance function parameters by cross
 validation and also by restricted maximum likelihood. For Kriging
 there is an easy to use function that also estimates the correlation
 scale (range parameter).  A major feature is that any covariance function
 implemented in R and following a simple format can be used for
 spatial prediction. There are also many useful functions for plotting
 and working with spatial data as images. This package also contains
 an implementation of sparse matrix methods for large spatial data
 sets and currently requires the sparse matrix (spam) package. Use
 help(fields) to get started and for an overview.  The fields source
 code is deliberately commented and provides useful explanations of
 numerical details as a companion to the manual pages. The commented
 source code can be viewed by expanding  source code version
 and looking in the R subdirectory. The reference for fields can be generated
 by the citation function in R and has DOI <doi:10.5065/D6W957CT>. 
  Development
 of this package was supported in part by the National Science Foundation  Grant
 1417857 and the National Center for Atmospheric Research. 
 
 
 Please see the excellent  vignette by Ashton Wiens and Mitch Krock on
 this package.
 [fieldsVignette](fieldsVignette.pdf)
 
#Browsing source code and help files
See the directory  [currentHelp](currentHelp) to browse the current versions of **fields** help files in html format. 

See the  directory [fields/R](fields/R) to browse the source code in
this package.
Note that orignal source code has numerous comments.
Unfortunately these are stripped out in the CRAN distribution and standard installation. 

#About the DOI version

The doi:10.5065/D6W957CT for the fields package is linked to the specific package version 8.4-1 

[fields_8.4-1.tar.gz](DOIinfo/fields_8.4-1.tar.gz)

MD5 check sum: 
 c89a47b9717d93ca9470a104bc8de088 
  This is the R format [DESCRIPTION](DOIinfo/DESCRIPTION) file distributed with this package 8.4-1.
#  
 Although this DOI version has the overall set of functions that are core to fields we caution the user that subsequent versions have fixed bugs and also added features that make the functions easier to use and more efficient. 
 

For the most recent version of fields please use a CRAN site such as the  [R studio CRAN  mirror site](http://cran.rstudio.com/) to download and install fields. Use citation("fields") in R to generate a citation with the current version number and this doi.





