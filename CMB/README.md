# NucTools Cluster Maps Builder: aggregate profile and aligned occupancy matrix visualizer

## Disclaimer

Cluster Maps Builder (CMB) is a part of NucTools 2.0 release. It is
at constant development and therefore may have some instability and bugs. Some known issues discribed in attached user manual.
If you find new bugs please contact y.vainshtein@zmbh.uni-heidelberg.de

## Introduction

Cluster Maps Builder (CMB) is primarily designed to visualize nucleosome
occupancy profiles of thousands of features aligned at genomic coordinate
corresponding to a specific feature, like transcription factor binding site or
transcription/translation initiation or termination site, using a heatmap
representation. CMB includes a K-means clustering step and is able to apply the
sorting/clustering order from initial matrix to a different matrix of the same size
and dimensions.

CMB is written using MATLAB and is using Java-based GUI (GUIDE). It can be downloaded as an executable file compiled to run under Windows operating system or as a sorce files.
The initial development was done using MATLAB 2014b but the program was
tested for compatibility with 2015a/b and 2016a. CMB has been tested for
Windows and MacOS X. It has not been tested for Linux yet.
