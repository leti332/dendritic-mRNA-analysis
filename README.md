# Dendritic mRNA analysis
Purpose: Analyze messenger RNA kinetics in neuronal dendrites

## Overview of Analysis Pipeline
* I. Image Pre-processing
* II. mRNA/Single Particle Tracking (SPT) Tracking
* III. SPT Post-processing and Analysis

For a detailed execution of this analytical framework, please see the SPTanalysis.md

## I. Image Pre-processing
<ins>Overview</ins>
* A. Create Average Z-projection using Kymograph Analysis Scripts (MATLAB)
* B. Generate dendritic outlines and crop movie (FIJI)

### A. Generate Average Z-projection using Kymograph Analysis Script (MATLAB)
Software Requirements:
* MATLAB 2018b
* Files
  - Kymograph_Analysis.m
  - kymo.m
  - avg_proj.m

Overview of Scripts
* compress 3D images into 2D images
* it generates an average of the compressed dimension
* the compressed dimension is y, since dendrites have limited height due to the narrow diameter

Directions
1. Have all files in the current path
2. Start with Kymograph_Analysis.m
  * Adjust input and output folders as needed
  * Run this file

### B. Generate dendritic outlines and crop movie (FIJI)
Software Requirements
* FIJI
  * ImageJ
* Files (ImageJ Macros)
  * Dendritic_Outlines.ijm
  * Segment_Dendrites_Movies.ijm
  * CombineDendrites.ijm
  * RollingBall-5px_Filter.ijm

Overview of Scripts
* *Dendritic_Outlines.ijm* - Saves outline of segment
* *Segment_Dendrites_Movies.ijm* - crops outline from movie
* *CombineDendrites.ijm* - concatenates cropped movies together
* *RollingBall-5px_Filter.ijm* - applies rolling ball filter with a diameter of 5 pixels

*Note: Future simplification possible.
there is an easier way to save outlines in the ROI manager of FIJI.
I realized this after I generated all these scripts.
It's likely this could be simplified by combining the ROI manager with the current scripts.*

Directions
* Open text editor in FIJI. Run scripts in the following order
1. Dendritic_Outlines.ijm
2. Segment_Dendrites_Movies.ijm
3. CombineDendrites.ijm
4. RollingBall-5px_Filter.ijm

*Note: The bandpass filter was selected since it yielded the best tracking results.
However, other filters for other datasets may be superior. For example, I also tried
bandpass filters, but got less accurate tracking.*

## II. Single Particle Tracking
<ins>Overview</ins>
* A. Diatrack - compatible with MSD

### A. Diatrack
Overview
* see this review article for directions
  * Vallotton P, van Oijen AM, Whitchurch CB, et al. Diatrack particle tracking software: Review of applications and performance evaluation. Traffic. 2017;18(12):840-852. doi:10.1111/tra.12530

Software Requirements
* MATLAB

Saving files
* save as ".MAT". This is needed for downstream processing in R


## III. SPT Post-processing and Analysis
### A. Diatrack SPT post-processing and analysis
#### i. R post-processing using modified scripts from [*sojourner* package ](https://github.com/sheng-liu/sojourner)
#### ii. MSD Analysis using modified scripts from [Lerner et al., 2020.](https://data.mendeley.com/datasets/hxnhtttxpk/1)
