# Dendritic mRNA analysis
Purpose: Analyze messenger RNA kinetics in neuronal dendrites

## Overview of Analysis Pipeline
* I. Image Pre-processing
* II. Single Particle Tracking
* III. SPT Post-processing and Analysis


## I. Image Pre-processing
<ins>Overview</ins>
* A. Create Average Z-projection using Kymograph Analysis Scripts (MATLAB)
* B. Generate dendritic outlines and crop movie (FIJI)
* C. Filter Movies to enhance mRNA detection

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


**Kymograph_Analysis.m**
```MATLAB
%% Clear All Variables before Beginning
clear all % clear all variables from workspace
close all % close all figures
clc % clear command window/console
%% Select Directory (Folder)
rootFolder =  'C:\\Users\\letinunez\\Desktop\\'; % directs GUI to this folder first
myDir = uigetdir(rootFolder, 'Select Folder with Images for Kymograph Analysis'); % uses GUI to get directory
% myDir = '/Users/letinunez/Desktop/ANALYSIS/KO/'; %MAC
% myDir ='C:\Users\lnunez\Desktop\ANALYSIS\KO\'; %PC - Lab
% myDir = 'C:\\Users\\letinunez\\Desktop\\Day2_KO20140717_WT20140716'; %PC- HOME
myFiles = dir(fullfile(myDir,'*.tif')); %gets all txt files in struct
%% Save Files
% Output Location for Saved Files
OutputPath = uigetdir(rootFolder, 'Select Output Folder for Kymographs'); % uses GUI to get directory
% OutputPath = '/Users/letinunez/Desktop/ANALYSIS/kymograph/'; %MAC
% OutputPath ='C:\Users\lnunez\Desktop\ANALYSIS\kymographs\'; % PC - Lab
% OutputPath = 'C:\\Users\\letinunez\\Desktop\\Kymographs'; %PC - HOME
%% Loop through Files in Directory (Folder)
% Progress Bar for Loop
f = waitbar(0, 'Starting');
n = length(myFiles);
% Loop through Files in myDir
for k = 1:n
  baseFileName = myFiles(k).name; %selects file from folder
  fullFileName = fullfile(myDir, baseFileName); %gets full path of file
  %fprintf(1, 'Now reading %s\n', fullFileName);
  % Runs Kymo function to generate Kymograph
  kymo(myDir, baseFileName, OutputPath)
  % Progress Bar
  waitbar(k/n, f, sprintf('Progress: %d %%', floor(k/n*100)));
  pause(0.1);
end
close(f) %closes progress bar
fprintf('Finished \n')
```
**kymo.m**
```MATLAB
function kymo(InputPath, filename, OutputPath)
%%
% fprintf(['InputPath: ' InputPath '\n'])
% fprintf(['File name: ' filename '\n'])
% Concatenate InputPath and filename to give imread full path of file
%pathfileName = [InputPath,filename]; % MAC
pathfileName = [InputPath,'\\', filename];% PC - HOME
% Read Image into MATLAB
testImage = imread(pathfileName,'tif',1);
% Print in Console Update
fprintf(['Reading Image: ' filename '\n']);
%% Get dimensions and stack number of movie
% Gets dimensions of the Image (read into Matlab)
[size1, size2] =size(testImage); %size1 - x pixels and size2 - y pixels
% Use imfinfo to get the number of stacks in the image
tiff_info = imfinfo(pathfileName);
TstackNum=length(tiff_info);
%fprintf(['Number of stacks: %f' TstackNum '\n']); % doesn't print; should
%try sprintf()
%% Create Average Projection
fprintf('Creating Average Projection\n')
avg_proj(TstackNum,OutputPath,pathfileName,1);
%% Create Kymograph
fprintf('Creating Kymograph from Image Stack\n')
Kymo=zeros(TstackNum,size2);
for i =1:TstackNum
    testImage=imread(pathfileName,'tif',i);
    Kymo(i,:)=squeeze(max(testImage,[],1));
end
%% Creat Kymograph Figure
fprintf('Creating Kymograph Figure\n')
figure('Name', 'Kymograph Analysis');
imagesc(Kymo)
Kymo=uint16(Kymo);
baseFileName = sscanf(filename,'%c',size(filename,2)-4); % may need to adjust
xlabel('X (pixels)')
ylabel('Time (frames)')
title([baseFileName ' ' 'Kymograph'])
%% Save Kymograph
fprintf(['Saving ' baseFileName '\n'])
% Save kymograph (_kymo,tif)
%imwrite(Kymo, [OutputPath, baseFileName, '_kymo.tif'], 'tif', 'WriteMode', 'overwrite', 'Compression', 'none'); % MAC
imwrite(Kymo, [OutputPath,'\\', baseFileName, '_kymo.tif'], 'tif', 'WriteMode', 'overwrite', 'Compression', 'none'); % PC
% Save kymograph with t and x units (_kymo_units.tif)
%saveas(gcf, [OutputPath, baseFileName, '_kymo_units.tif'])%MAC
saveas(gcf, [OutputPath,'\\', baseFileName, '_kymo_units.tif']) %PC
% Save matlab figure (_kymo.fig)- NOT WORKING
%saveas(gcf, [OutputPath, baseFileName, '_kymo.fig'], '.fig')
end
```
**avg_proj.m**
```Matlab
%% Creates Average Z-projection of movie
function OutputImage=avg_proj(ZstackNum,output_PathName,PathFileName,option)  % StartImageNum >10,  EndImageNum<100
%% Define Path and File Name
[pathName, fileName] = fileparts(PathFileName); % get path and filename (without extension) separately
%PathName=[pathName,'\']; %MAC
PathName = [pathName, '\\']; % PC
FileName=[fileName,'.tif'];
%% Read image into MATLAB
testImage = imread(PathFileName);
%% Generate Average Projection
% Allocate memory
ch1 = uint16(zeros(size(testImage,1),size(testImage,2),ZstackNum));
proj = uint16(zeros(size(testImage,1),size(testImage,2)));
    for i =1:ZstackNum
        ch1(:,:,i) = imread(PathFileName, i);
    end
proj(:,:) = mean(ch1,3);
%ch1_proj = maxproj(ch1);%ch1_proj=uint16(ch1_proj);
OutputImage=proj;
%%
if option == 1
%    imwrite(proj, [output_PathName,'AVG_',fileName, '.tif'], 'tif', 'WriteMode', 'append', 'Compression', 'none');%MAC
    imwrite(proj, [output_PathName,'\\','AVG_',fileName, '.tif'], 'tif', 'WriteMode', 'append', 'Compression', 'none');%PC
end
```
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

*Note: The rolling ball filter was selected since it yielded the best tracking results.
However, other filters for other datasets may be superior. For example, I also tried
bandpass filters, but got less accurate tracking.*

**Dendritic_Outlines.ijm**
```javascript
// Save Dendritic Outlines
// WORKS - 18MAR2021
// Editor: Leti Nunez
//////////////////////////////////////////////////////////////////

// 1. Load Average Z-projection (see Kymograph Analysis MATLAB script) and raw movie
// 2. Use Segmented Line on Avg Z-project to segment dendrites
// Note: Start from soma (if no soma in movie, base on branching direction. Bifurcating
// branching usually are in the anterograde direction. In other words, the "Y" should
// be to the right when segmenting.
// 3. Run this script to save dendritic outline


// Specify the dendrite number
	dNum = "-" + "2";

// Apply dendritic outline to duplicated image (for saving purposes)
	orig = getImageID(); //get imageid of selected image
	//print(title);//sanity check
	run("Duplicate...", " "); // make duplicate image it can be saved with outline

	new = getImageID(); // save imageid of duplicate image

	selectImage(orig); //select original image to get the outline selection
	selectImage(new); //select duplicate image
	run("Restore Selection"); // apply outline to duplicate image

// Change Title to reflect dendrite number
	selectImage(orig); //select original image for title
	title = getTitle();
	idx = indexOf(title, "."); // get index of period
	//print(idx);//sanity check

	Title = substring(title, 4, idx); //get title before the period
	print(Title);//sanity check


	Title = "Outline_" + Title + dNum + ".tif"// save outline
	//print(Title);//sanity check

// Get Output Directory for Image
//out = getDirectory("Output Path"); // GUI selects directory
	out = "C:/Users/letinunez/Dropbox (EinsteinMed)/Projects/RNA-Protein_Dynamic_Modeling/ACTBmRNAtracking_ZBP1KO/PAPER_FINAL/Data/2_Dendrite Outlines/WT/exp3_20150717"
	//print(out);//sanity check

// Save Image
	selectImage(new);
	saveAs("Tiff", out + File.separator + Title);
```
**Segment_Dendrites_Movies.ijm**
```javascript
#@ File(label="Raw Movie directory", style = "directory") movies
#@ File(label="Dendritic Outline directory", style = "directory") outlines
#@ File(label="Segmented Dendrite Output directory", style = "directory") output
#@ String(label = "File Suffix", value = ".tif") suffix


processFolder(movies, outlines);

function processFolder(movies, outlines) {
	print("Processing...");
	mlist = getFileList(movies);

	olist = getFileList(outlines);

	for (i = 0; i < mlist.length; i++) {
		if(File.isDirectory(movies + mlist[i]))
			processFolder("" + movies + mlist[i]);
			//print("Movie: " + mlist[i]);
			idx = indexOf(mlist[i], suffix); //get index of suffix
			Movie = substring(mlist[i],0, idx);// get name without suffix
			//print("Movie: " + Movie); //print name without suffix
		if(endsWith(mlist[i], suffix)){
			//Loop through outline folder
			for (k = 0; k < olist.length; k++){
				if(File.isDirectory(outlines + olist[k]))
					processFolder("" + outlines + olist[k]);
					//print("Outline: " + olist[k]);

				// Conditional if the Movie name is a substring of the outline name
				//Execute following code
				if(indexOf(olist[k], Movie) >= 0){
					strDend(outlines, movies, olist[k], mlist[i], output);
				}else{
					print(olist[k] + " & " + mlist[i] + " are not Matches");
				}
			}
		}
	}
	print("Finished!");

}

function strDend(outlines, movies, outline, movie, output){
	setBatchMode(true); // true prevents image windows from opening while the script is running

	print("Matches: " + outline + " & " + movie);
	open(outlines + File.separator + outline);
	print("Processing: " + outlines + File.separator + outline);
	oId = getImageID();
	otitle = getTitle();


	open(movies + File.separator + movie);
	mId = getImageID();
	mtitle = getTitle();
	selectImage(oId);
	selectImage(mId);
	run("Restore Selection");
	run("Straighten...", "title line=50 process");

	// Remove string "Outline_" from title name
	rm = "Outline_";
	index = lengthOf(rm);

	ix = lastIndexOf(otitle, "Outline_");
	index = index + ix;

	file = substring(otitle,index);

	// Save File
	print("Saving Dendritic: " + file);
	saveAs("tiff", output + File.separator + file);
}
```

**CombineDendrites.ijm**
```javascript
#@ File(label="Dendrite directory", style = "directory") movies
#@ File(label="Combined Dendrite Output directory", style = "directory") output
#@ String(label = "File Suffix", value = ".tif") suffix


processFolder(movies);

function processFolder(movies) {
	print("Processing...");
	mlist = getFileList(movies);

	for (i = 0; i < mlist.length; i++) {

		file = "KO_Comb_20140717_KO23.tif";

		if(endsWith(mlist[i], suffix)){
			if(i == 0){
				print(i);
				comb1stDend(movies, mlist[i], mlist[1], output, file);
			}
			if(i > 1){
				print("Progress: " + i + "/" + mlist.length);
				combDend(movies, file, mlist[i], output);
			}

		}
	}
	print("Finished!");
}

function comb1stDend(movies, movie1, movie2, output, file){
	setBatchMode(true); // true prevents image windows from opening while the script is running

	//Movie1
	open(movies + File.separator + movie1); //open file
	print("Movie1: " + movie1);
	idx = indexOf(movie1, suffix); //get index of suffix
	Movie1 = substring(movie1,0, idx);// get name without suffix
	//M1title = getTitle();
	//print("Movie1: " + M1title);

	//Movie2
	open(movies + File.separator + movie2); //open file
	print("Movie2: " + movie2);
	idx = indexOf(movie2, suffix); //get index of suffix
	Movie2 = substring(movie2,0, idx);// get name without suffix

	//Combine Dendrites
	run("Combine...", "stack1=" + movie1 + " stack2=" + movie2 + " combine");

	print(output);
	saveAs("tiff", output + File.separator + file);

	close();
	print("1st Combination Completed");
}

function combDend(movies, movie1, movie2, output){
	setBatchMode(true); // true prevents image windows from opening while the script is running

	//Movie1 - Combined Movie
	open(output + File.separator + movie1); //open file
	print("Movie1: " + movie1);
	idx = indexOf(movie1, suffix); //get index of suffix
	Movie1 = substring(movie1,0, idx);// get name without suffix
	//M1title = getTitle();
	//print("Movie1: " + M1title);

	//Movie2
	open(movies + File.separator + movie2); //open file
	print("Movie2: " + movie2);
	idx = indexOf(movie2, suffix); //get index of suffix
	Movie2 = substring(movie2,0, idx);// get name without suffix

	//Combine Dendrites
	print("Combining: " + movie1 + " & " + movie2);
	run("Combine...", "stack1=[" + movie1 +"] stack2="+ movie2 +" combine");

	saveAs("tiff", output + File.separator + movie1);
	close();
	print("Completed");
}
```
### C. Filter Movies to enhance mRNA detection
Overview: To enhance mRNA detection, segmented dendrites are filtered. We recommend trying a number of filters first to empirically determine which is the ideal filter for your dataset. Below are two scripts that apply filters. The same filter should be applied to all data to avoid artefacts introduced by changing filters. 

Files (use one or the other)
- Process_Folder_BandpassFFT_3-5px.ijm
- RollingBall-5px_Filter.ijm

**Process_Folder_BandpassFFT_3-5px.ijm**
``js
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix


processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	setBatchMode(true); // true prevents image windows from opening while the script is running
	open(input + File.separator + file); //open file
	
	print("Processing: " + input + File.separator + file);
	run("Bandpass Filter...", "filter_large=5 filter_small=3 suppress=None tolerance=5 autoscale saturate process");
    run("Enhance Contrast", "saturated=0.35");
	
	print("Saving to: " + output);
	saveAs("tiff", output + File.separator + "bandpass5-3px_" + file);
}
``
**RollingBall-5px_Filter.ijm**
```js
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
// See also Process_Folder.py for a version of this code
// in the Python scripting language.
processFolder(input);
// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	print("Finished!");
}
function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	setBatchMode(true); // true prevents image windows from opening while the script is running
	open(input + File.separator + file); //open file

	print("Processing: " + file);
	//run("Bandpass Filter...", "filter_large=4 filter_small=2 suppress=None tolerance=5 autoscale saturate process");
    //run("Enhance Contrast", "saturated=0.35");
    rpx = "5";
	run("Subtract Background...", "rolling=" + rpx +" stack");

	print("Saving");
	//saveAs("tiff", output + File.separator + "bandpass4-2px_" + file);
	saveAs("tiff", output + File.separator + "rollingball" + rpx + "px_" + file);
}
```

## II. Single Particle Tracking
<ins>Overview</ins>
* A. Diatrack

### A. Diatrack
Overview
* see this review article for directions
  * Vallotton P, van Oijen AM, Whitchurch CB, et al. Diatrack particle tracking software: Review of applications and performance evaluation. Traffic. 2017;18(12):840-852. doi:10.1111/tra.12530

Software Requirements
* MATLAB

Saving files
* save as ".MAT". This is needed for downstream processing in R


## III. SPT Post-processing and Analysis
### A. mRNA track post-processing
Overview: convert Diatrack mRNA trajectories into a common format that can utilized for trajectory analyses
Software
* R
Files
- Download all dependencies for sojourner package
- DiatrackAnalyzerBYCELL_v2.Rmd
#### i. R post-processing using *sojourner* package

**DiatrackAnalyzerBYCELL_v2.Rmd**
###### 1. Import Data

* must have folder with all matlab files together

###### 2. Load sojourner package
```r
# load library
library(sojourner)
```
###### 3. Specify Data Type for Import
```r
# specify condition value
# WT = 0
# KO = 1
condition = 1
# conditional statement assigns strings to condition value
if(condition == 0){
  print('Processing WT Files')
  condition = "WT"
}else if(condition == 1){
  print('Processing ZBP1KO Files')
  condition = "ZBP1KO"
}
```
###### 4. Specify Folders for import
```r
# Specify Folders
#exp1 = dirname(file.choose())
exp1 ="/Users/letinunez/Dropbox (EinsteinMed)/Projects/RNA-Protein_Dynamic_Modeling/ACTBmRNAtracking_ZBP1KO/PAPER_FINAL/Data/6_Track/Diatrack/BYCELL/KO/20140503_Rollingball_250pxlowerlimit_max5pxjump"
#exp2 = dirname(file.choose())
exp2 = "/Users/letinunez/Dropbox (EinsteinMed)/Projects/RNA-Protein_Dynamic_Modeling/ACTBmRNAtracking_ZBP1KO/PAPER_FINAL/Data/6_Track/Diatrack/BYCELL/KO/20140710_Rollingball_150pxlowerlimit_max5pxjump"
#exp3 = dirname(file.choose())
exp3 = "/Users/letinunez/Dropbox (EinsteinMed)/Projects/RNA-Protein_Dynamic_Modeling/ACTBmRNAtracking_ZBP1KO/PAPER_FINAL/Data/6_Track/Diatrack/BYCELL/KO/20140716_Rollingball_150pxlowerlimit_max5pxjump"
# print Processing Diatrack Files
print('Processing Diatrack trajectories by CELL')
t1 = createTrackll(exp1, input=2, ab.track=FALSE, cores=4)
t2 = createTrackll(exp2, input=2, ab.track=FALSE, cores=4)
t3 = createTrackll(exp3, input=2, ab.track=FALSE, cores=4)
```
###### 5. Save trajectory info as RData file
```r
file = paste(condition,exp,"tracks.RData", sep ="-")
save(t1, file = file)
exp = "exp2"
file = paste(condition,exp,"tracks.RData", sep ="-")
save(t2, file = file)
exp = "exp3"
file = paste(condition,exp,"tracks.RData", sep ="-")
save(t3, file = file)
```
###### 6. Set track to look at
```r
exp = "exp1"
tracks = t1
```
###### 7. Look at Names of Cells
```r
# Features of Imported Data
# Condition Names = Name of Folders
names(tracks) #
# If needed, rename folders/conditions
# names(tracks) = c("exp1", "exp2", "exp3") # Names of Folders
# Cell Names
length(names(tracks))
```
###### 8. Filter Data
1. Track Length
- remove trajectories which are smaller than 15 consecutive frames and greater than 200 consecutive frames

2. Pixel Intensity
- remove trajectories which are smaller than 200 pixels and greater than 1000 pixels

```r
#--------------------------------------------------------------------
# Filter Tracks based on Length
## function
filterTrack = function(trackll, filter = c(min = 7, max = Inf)) {

    # reinforce name
    names(filter) = c("min", "max")

    cat("applying filter, min", filter["min"], "  max", filter["max"],
        "\n")

    track.len = list()
    for (i in seq_along(trackll)) {
        track.len[[i]] = vapply(trackll[[i]], function(track) {
            dim(track)[1]
        }, integer(1))
        trackll[[i]] = trackll[[i]][track.len[[i]] >=
            filter["min"] & track.len[[i]] < filter["max"]]
    }

    return(trackll)
}
#--------------------------------------------------------------------
#Filter Tracks based on Intensity
## function
filterInten = function(trackll, filter = c(min = 200, max = Inf)) {

    # reinforce name
    names(filter) = c("min", "max")

    cat("applying filter, min", filter["min"], "  max", filter["max"],
        "\n")

    track.maxint = list()
    for (i in seq_along(trackll)) {
        track.maxint[[i]] = lapply(trackll[[i]], function(track) {
            max(track[5])
        })
        trackll[[i]] = trackll[[i]][track.maxint[[i]] >=
            filter["min"] & track.maxint[[i]] < filter["max"]]
    }

    return(trackll)
}

#--------------------------------------------------------------------
# Apply Filters
trackll_f <- filterTrack(tracks)
trackll_i <- filterInten(tracks)
trackll_i_low = filterInten(tracks, c(100,375))
trackll_fi <- filterTrack(trackll_i)

file = paste(condition,exp,"trackfiltered.RData", sep ="-")
save(trackll_f, file = file)

# # KO only - Filtering
# trackll_KO_f <- filterTrack(trackll_KO) # 15 to 200
# trackll_KO_i <- filterInten(trackll_KO) #200 to 1000
# trackll_KO_fi <- filterTrack(trackll_KO_i, filter = c(min = 15, max = Inf))
```
###### 9. Calculate Maximum Intensity for Each Trajectory/Particle
```r
#--------------------------------------------------------------------
# calculate the maximum intensity for each particle - WOrks
maxIntensity = function(track){
  track.maxint = list()
  for (i in seq_along(track)){

  track.maxint[[i]] = lapply(track[[i]], function(x){
    max(x[5])
  })

  }
  names(track.maxint) = names(track)
return(track.maxint)
}
#--------------------------------------------------------------------
# Apply function to calculate maximum intensities
maxInten_tracks = maxIntensity(tracks) # unfiltered
maxInten_trackllf = maxIntensity(trackll_f) # Length filtered
maxInten_trackllfi = maxIntensity(trackll_fi) # Length & Intensity filtered

maxInten_tracklli_LOW = maxIntensity(trackll_i_low)
#save(maxInten_tracks, maxInten_trackllf, file = "maxInten_rb5px.RData")
```
###### 10. Reshape Intensity Data for plotting
use <r><reshape2::melt()> to wrangle data more easily
```r
# Reshape Data
maxInten_melt_unfilt = reshape2::melt(maxInten_tracks) # unfiltered data
maxInten_melt_f = reshape2::melt(maxInten_trackllf) # Length filtered

maxInten_melt_fi = reshape2::melt(maxInten_trackllfi) # Length & intensity filtered
maxInten_melt_i_low = reshape2::melt(maxInten_tracklli_LOW)
# Rename column names for plotting
colnames(maxInten_melt_unfilt) = c("value", "mRNP", "cell")
colnames(maxInten_melt_f) = c("value", "mRNP", "cell")
colnames(maxInten_melt_fi) = c("value", "mRNP", "cell")
colnames(maxInten_melt_i_low) = c("value", "mRNP", "cell") # dim particles
```
###### 11. Plot Maximum Intensities of Trajectories
```r
library(ggplot2)
intensity.plot = function(dataframe,title,x.scale = c(min=0,max=1000), y.scale = c(min=0,max=1000)) {
  # Histogram plot of Intensities
  ggplot(dataframe, aes_string(
    x = "value"
    ,group = "cell"
    ,fill = "cell"
    )) +
  geom_histogram(position = "dodge"
                 #, colour = "white"
                 ) +
  xlim(x.scale["min"], x.scale["max"])  +
  ylim(y.scale["min"], y.scale["max"]) +
  labs(
    x = "Maximum Intensity of mRNPs",
    y = "Number of mRNPs") +
  ggtitle(title)

}
intensity.plot(maxInten_melt_unfilt,"mRNPs Tracked from Diatrack", x.scale = c(min=200,max=1000), y.scale = c(min=0,max=1000))
intensity.plot_f = intensity.plot(maxInten_melt_f,"Length Filtered Data")
intensity.plot_fi = intensity.plot(maxInten_melt_fi,"Length & Intensity Filtered Data")
intensity.plot_i_low = intensity.plot(maxInten_melt_i_low,"Dim mRNP")
intensity.plot_unfilt
intensity.plot_f
intensity.plot_fi
intensity.plot_i_low
# Raw Data Summary
print("Raw Data")
summary(maxInten_melt_unfilt)
summary(maxInten_melt_unfilt[maxInten_melt_unfilt$condition=="WT",])
summary(maxInten_melt_unfilt[maxInten_melt_unfilt$condition=="ZBP1KO",])
# Filtered Summary
print("Length Filtered Data")
summary(maxInten_melt_f)
print("Length & Intensity Filtered Data")
summary(maxInten_melt_fi)
summary(maxInten_melt_fi[maxInten_melt_fi$condition=="WT",])
summary(maxInten_melt_fi[maxInten_melt_fi$condition=="ZBP1KO",])
print("Dim mRNPs")
summary(maxInten_melt_i_low)
summary(maxInten_melt_i_low[maxInten_melt_i_low$condition=="WT",])
summary(maxInten_melt_i_low[maxInten_melt_i_low$condition=="ZBP1KO",])
#--------------------------------------------------------------------
#t.interval = 100 # time interval in milliseconds
t.interval = 0.100 # time interval in seconds
#--------------------------------------------------------------------
# function to calculate time values
.dwellTime = function(trackl, t.interval = 10) {
    vapply(trackl, function(x) {
        dim(x)[1] * t.interval
    }, FUN.VALUE=double(1))
}
#--------------------------------------------------------------------
# Apply function using lappy to get dwelltime for each trajectory
dwell.time_unfilt = lapply(tracks, function(x) {
        .dwellTime(x, t.interval)
    })
dwell.time_f = lapply(trackll_f, function(x) {
        .dwellTime(x, t.interval)
    })

dwell.time_fi = lapply(trackll_fi, function(x) {
        .dwellTime(x, t.interval)
    })
dwell.time_i_low = lapply(trackll_i_low, function(x) {
        .dwellTime(x, t.interval)
    })
# Reshape Data for plotting
dwell.time.mlt_unfilt = reshape2::melt(dwell.time_unfilt)
dwell.time.mlt_f = reshape2::melt(dwell.time_f)
dwell.time.mlt_fi = reshape2::melt(dwell.time_fi)
dwell.time.mlt_i_low = reshape2::melt(dwell.time_i_low)
save(dwell.time.mlt_unfilt, dwell.time.mlt_f, file = "dwellTime_bp3-5px.RData")
save(dwell.time.mlt_i_low, file = "dwellTime_rb5px_dim.RData")
# Rename columns for plotting
colnames(dwell.time.mlt_unfilt) = c("value", "variable")
colnames(dwell.time.mlt_f) = c("value", "variable")
colnames(dwell.time.mlt_fi) = c("value", "variable")
colnames(dwell.time.mlt_i_low) = c("value", "variable")
```
###### 12. Plot Data
```r
dwellTime.plot = function(dataframe, title, x.scale = c(min = 0,max = 2.5), y.scale = c(min = 0,max = 1000)){
  ggplot(dataframe, aes_string(x = "value",
        group = "variable", fill = "variable")) +
        geom_histogram(binwidth = t.interval, position = "dodge") +
    xlim(x.scale["min"], x.scale["max"]) +
    theme_bw() +
    ylim(y.scale["min"], y.scale["max"]) +
        theme(legend.title = element_blank()) +
        labs(x = "Lifetime of mRNP tracks (sec)", y = "Number of mRNPs") +
ggtitle(title)
}
dwellTime.plot(dwell.time.mlt_unfilt, "Raw Data",y.scale = c(min=0,max=250))
dwellTime.plot(dwell.time.mlt_f, "Length Filtered Data", x.scale = c(min=1.3,max=2.5), y.scale = c(min=0,max=1000))
dwellTime.plot(dwell.time.mlt_fi, "Length & Intensity Filtered Data", x.scale = c(min=1.3,max=2.5))
dwellTime.plot(dwell.time.mlt_i_low, "Dim mRNPs")
```
###### 13. Summarize Trajectory Lifetimes
- Compare Filtered and unfiltered/Raw data

```r
# Raw Data Summary
print("Raw Data")
summary(dwell.time.mlt_unfilt)
# Filtered Summary
print("Length Filtered Data")
summary(dwell.time.mlt_f)
# Filtered Summary
print("Length & Intensity Filtered Data")
summary(dwell.time.mlt_fi)
save(dwell.time_unfilt, dwell.time_f, dwell.time_fi, file = "dwelltime.RData")
```
###### 14. Export Trajectories using built-in sojourner functions
```r
exportTrackll(tracks, cores = 4)
exportTrackll(trackll_f, cores = 4)
exportTrackll(trackll_fi, cores = 4)
exportTrackll(trackll_i_low, cores = 4)
```
###### 15. Read in CSVs and export into TPM csv format
```{r}
#myDir = dirname(file.choose())
myDir ="/Users/letinunez/Dropbox (EinsteinMed)/Projects/RNA-Protein_Dynamic_Modeling/ACTBmRNAtracking_ZBP1KO/PAPER_FINAL/Scripts/CSV"
myFiles = list.files(myDir)
num_files = length(myFiles)
#OutputPath = dirname(file.choose())
OutputPath = "/Users/letinunez/Dropbox (EinsteinMed)/Projects/RNA-Protein_Dynamic_Modeling/ACTBmRNAtracking_ZBP1KO/PAPER_FINAL/Data"
t.interval = 0.1
# Convert Frames to Seconds
convertFrames2sec = function(listOftraj){
  newTracks = listOftraj
  for (i in seq_along(listOftraj)){
  newTracks[[i]]$t = vapply(listOftraj[[i]]$Frame,function(x){
     as.numeric(x*t.interval)
   }, FUN.VALUE = 1.0)

  }
  names(newTracks) = names(listOftraj)
return(newTracks )
}
## Make List with all Trajectory Information from Folder
df = list()
for (i in seq_along(myFiles)) {
  print(myFiles[i])
  file = paste(myDir,myFiles[i], sep = "/")
  df[[i]] = read.csv(file)
}
names(df) = myFiles
df_sec = convertFrames2sec(df)
head(df_sec[[1]])
vars = c("Trajectory", "Frame", "t", "x", "y")
## Run through each file separately
for (i in seq_along(df_sec)){
 # print(i)
    TPM = dplyr::select(df_sec[[i]], one_of(vars))
    colnames(TPM) = c("trajectory", "frame", "t", "x", "y")
    col_order = c("frame", "t", "trajectory", "x", "y")
    data_TPM = TPM[,col_order]
    fileName = strsplit(names(df_sec[i]), ".csv")
    fileName = paste(fileName, "TPM.csv", sep = "-")
    fullFile = paste(OutputPath,fileName, sep ="/")
    #print(fullFile)
    write.csv(data_TPM, fullFile)
}
```
### B. mRNA track analysis using MSD based on Lerner et al. 2020
Overview: Stratefy different motion types (e.g. confined vs. directed) using mean square displacement analysis based on Lerner et al. 2020. 
Software
* MATLAB
Files
- Download all scripts from Lerner et al. 2020
- skeleton.m
- ClassifyMotion.m
- MSD_LinearFit.m
	MSDAnalysis.m
- MSD_ParabolicFit.m
Directions
Run skeleton.m
Select MSD_LinearFit.m or MSD_ParabolicFit.m based on particle motion. 
	
##### 1. Run skeleton.m
**skeleton.m**
```Matlab
%%
%% Select Directory (Folder) with Trajectory Files
rootFolder =  'C:\\Users\\letinunez\\Desktop\\'; % directs GUI to this folder first
% myDir = uigetdir(rootFolder, 'Select Folder with trajectory files for MSD Analysis'); % uses GUI to get directory
% myDir = '/Users/letinunez/Desktop/ANALYSIS/KO/'; %MAC
% myDir = '/Users/letinunez/ACTBmRNAZBP1KOneurons/inst/Rmarkdown files/TPM_files'; %MAC
% myDir ='C:\Users\letinunez\Dropbox (EinsteinMed)\Projects\Project1-ZBP1KO\ACTBmRNAZBP1KOneurons\inst\Analysis\TPM_matlab_input\WT\'; %PC - Lab
myDir = 'E:\New Dropbox\Dropbox (EinsteinMed)\ACTBmRNAZBP1KOneurons\inst\Rmarkdown files\TPM_files';% PC- HOME
myFiles = dir(fullfile(myDir,'*.csv')); %gets all csv files in struct

%% Save Files

% Output Location for Saved Files
% OutputPath = uigetdir(rootFolder, 'Select Output Folder for MSD Results'); % uses GUI to get directory
% OutputPath = '/Users/letinunez/Desktop/ANALYSIS/kymograph/'; %MAC
% OutputPath = '/Users/letinunez/ACTBmRNAZBP1KOneurons/inst/Analysis_28APR2022'; %MAC
% OutputPath ='C:\Users\letinunez\Dropbox (EinsteinMed)\Projects\Project1-ZBP1KO\ACTBmRNAZBP1KOneurons\inst\Analysis\TPM_matlab_output\'; % PC - Lab
OutputPath = 'E:\New Dropbox\Dropbox (EinsteinMed)\ACTBmRNAZBP1KOneurons\inst\Analysis_29APR2022';% PC- HOME

%% Input Variables for MSD


%% Loop through Files in Directory (Folder)
% Progress Bar for Loop
f = waitbar(0, 'Starting');

n = length(myFiles); % number of files in folder

% Loop through Files in myDir
for k = 1:n
  baseFileName = myFiles(k).name; %selects file from folder
  fullFileName = fullfile(myDir, baseFileName); %gets full path of file
  fprintf(1, 'Now reading %s\n', fullFileName);



  % Run MotionClassifier
  fprintf('\nStarting Motion Classifier\n')
  ClassifyMotion(fullFileName,inputVar,OutputPath)
  close all
  % Progress Bar
  waitbar(k/n, f, sprintf('Progress: %d %%', floor(k/n*100)));
  pause(0.1);


end
close(f) %closes progress bar

fprintf('Finished \n')
```
**ClassifyMotion.m**
```Matlab
function ClassifyMotion(fullFileName,inputVar,OutputPath)

%% Clear all variables Except the following
clearvars -except baseFileName fullFileName f inputVar k myDir myFiles n OutputPath rootFolder
[pathName, baseName] = fileparts(fullFileName); % get path and filename (without extension) separately
%fprintf(strcat('Input Path: ', pathName))
fprintf(strcat('Processing : ',baseName,'\n'))

% Make Folder in Output Folder with Base Name
OutputFolder = fullfile(OutputPath, 'MSD_Results');
mkdir(OutputFolder);

%% Initialize Variables
% FILE TYPE
if inputVar.fileType == 0
   fprintf('\nProcessing TrackMate File...\n')
   TrackMate = 1; %If the track file is an .xml file from TrackMate.
   SlimFast = 0;
elseif inputVar.fileType == 1
   fprintf('\nProcessing SlimFast File...\n')
   TrackMate = 0;
   SlimFast = 1; %If the track file is a .csv file from slimFAST.
end

fprintf('File Parameters ------------------\n')
% SPACE UNITS
if inputVar.spaceunits == 0
   fprintf('Space units = pixels\n')
   Input_space_units = 'pixels';
elseif inputVar.spaceunits == 1
   fprintf('Space units = microns\n')
   Input_space_units = 'microns';
end

Frame_interval = inputVar.frameinterval;

% TIME UNITS
if inputVar.timeunits == 0
   fprintf('Time units = frames\n')
   Input_time_units = 'frames';
elseif inputVar.timeunits == 1
   fprintf('Time units = seconds\n')
   Input_time_units = 'seconds';
end

pixel_size = inputVar.pixelsize;

% TRAJECTORY DIMENSIONS
dimension = inputVar.dim;
if dimension == 2
   fprintf('Trajectories are 2D\n')
elseif dimension == 3
   fprintf('Trajectories are 3D\n')
end

% FILTER LENGTH
minimum_track_length = inputVar.FilterMinLen;
% fprintf('All trajectories less than %i will be filtered\n', minimum_track_length)

% CONFINED TRAJECTORY THRESHOLD
Jump_confined_threshold = inputVar.jumpConf;
% fprintf('\nConfined Trajectories will be separated into two groups by %.4f um [Average Jump]',Jump_confined_threshold)
% MSD FITTING PARAMETERS
% Linear Fitting
TMSD_fitting_points = inputVar.TMSDfitpt;
% fprintf('\nLinear Regression of T-MSD will use %i to fit Diffusion Coefficient with Confidence Intervals', TMSD_fitting_points)
if TMSD_fitting_points < 3
   fprintf('Linear Fitting Requires a minimum of 3 points for Confidence Interval calculations\n')
end
TEMSD_fitting_points = inputVar.TEMSDfitpt;  %Number of points used of the TE-MSD for fitting the Diffussion Coefficient. (LINEAR FITTING) The minimum 3 points for being able to calculate the Confidence Intervals
% fprintf('Linear Regression of T-MSD will use %i to fit Diffusion Coefficient with Confidence Intervals', TEMSD_fitting_points)
if TEMSD_fitting_points < 3
   fprintf('Linear Fitting Requires a minimum of 3 points for Confidence Interval calculations\n')
end
% Logarithmic Fitting
TLOGLOG_fitting_points = inputVar.Tloglogfitpt; %Number of points used of the TE-MSD for fitting the Diffussion Coefficient. (POWER-LAW FITTING).
% fprintf('\nLogarithmic Regression of T-LogLogMSD will use %i to fit Diffusion Coefficient with Confidence Intervals', TLOGLOG_fitting_points)

max_alpha_conf = inputVar.maxAlphaconf; %Maximum alpha value for the power-fitting of the TE-MSD in order to consider Confined Motion
% fprintf('\nPower-fitting of TE-MSD will define Confined Motion below an alpha value of : ', max_alpha_conf)
min_alpha_directed = inputVar.minAlphadir; %Minimum alpha value for the power-fitting of the TE-MSD in order to consider Directed Motion
% fprintf('\nPower-fitting of TE-MSD will define Directed Motion above an alpha value of : ', min_alpha_directed)
% fprintf('\nPower-fitting of TE-MSD will define Brownian Motion above an alpha value of %.4f and below an alpha value of %.4f', max_alpha_conf, min_alpha_directed)

% R-squared Threshold
R2LIMIT = inputVar.R2LIMIT;
% fprintf('\nAll Tracks with an R-squared value between 1 and %.4f will be used for downstream analysis', R2LIMIT)
%------------------------------------------------------------------
%Butterfly Trajectories---------------------------------------------------
%------------------------------------------------------------------
% fprintf('\nButterfly Trajectories Parameters----')
jump_threshold = inputVar.butterfly_jumpthreshold;
% fprintf('\nOne Track MUST jump more than its own (average_jump + %.4f *std_jump) to be considered a Butterfly track', jump_threshold)
minim_dist = inputVar.butterfly_minDist;
% fprintf('\nOne Track MUST travell a total distance bigger than its own average_jump*%.4f to be considered a Butterfly track', minim_dist)
Conf2JumpThr = inputVar.butterfly_Conf2JumpThr;
% fprintf('\nButterfly track need to have a Jump equal or bigger than the average radius of confinement of its confined segments multiplied by %.4f',Conf2JumpThr)
Out_percentage = inputVar.butterfly_OutPercentage;
% fprintf('\n%i Minimum Percentage of points that a JUMP must have OUTSIDE the previous and posterior polygon (CONVEXHULL) to be considered an OUTER Segment. (Number from 0 to 100).', Out_percentage)
N_sliding = inputVar.butterfly_Nsliding;
% fprintf('\n%i points to check linearity of segments',N_sliding)
P_min_linear = inputVar.butterfly_P_min_linear;
% fprintf('\nMinimum to consider that a segment is linear : %.4f',P_min_linear)
angleTH = inputVar.butterfly_angleTH;
% fprintf('\nMinimum Angle for considering a jump as "directed" in butterfly tracks (degrees) Direct(>angleTH) or non-direct(<angleTH) : %i',angleTH)
Max_Jump_Confined_Butterfly = inputVar.butterfly_MaxJumpConf;
% fprintf('\nMaximum jump a confined segment of a butterfly track can have. (Use for track segmentation). : %.4f um',Max_Jump_Confined_Butterfly)
Min_num_points_butt_confined_segment = inputVar.butterfly_MinNumPtConfSeg;
% fprintf('\n Minimum number of points that a confined segment of a butterfly track can have : %i. (The tracks that do not fulfill this condition will be discarded)', Min_num_points_butt_confined_segment)
%------------------------------------------------------------------
%Circle Confined Diffusion---------------------------------------------------
%------------------------------------------------------------------
% fprintf('\nCircle of Confined Diffusion Parameters----')
Offset0 = inputVar.ConfDiffOffset;
% fprintf('\n%.4f is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in ï¿½m]. ', Offset0)
num_points = inputVar.ConfDiffNum_points;
% fprintf('\n%i points will be used to fit the confined diffussion circle model to the TE-MSD.', num_points)
D0 = inputVar.ConfDiffD0;
% fprintf('\nStarting value for the least squares fitting for the Diffusion Coefficient (ï¿½m^2/s) on the Confined Circle Diffusion Model : %.4f',D0)
R0 = inputVar.ConfDiffR0;
% fprintf('\n%Starting value for the least squares fitting for the Radius of Confinement (ï¿½m) on the Confined Circle Diffusion Model : %.4f',R0)

%% Initialize Variables----------------------------------------
%%----------------------------------------------------------
%--------------------------------------------------------------
SPACE_UNITS = 'µm'; %This is just for visualization purposes
TIME_UNITS = 's'; %This is just for visualization purposes
ma = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);  %Initialize MSD analyzer
ma_AllTracks = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_confined = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_confined_High_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_confined_Low_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_directed = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_brownian = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments_confined = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments_directed = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
fprintf('\n---All variables initialized---\n')


%% Load the Data:
%-------------------------------------------------------------------------
%Tracks-------------------------------------------------------------
%--------------------------------------------------------------------------

% Not working
% if TrackMate == 1;
% %Extract the trajectories and other related information from the results of the tracking
%     [trajectory{1}, metadata] = importTrackMateTracks(fullFileName),'clipZ',1);
%
% end


%Load data from SlimFast
if SlimFast == 1
   [trajectory{1}, tracklength{1}] = slimFast_ImportTracks(fullFileName);
end

fprintf('\nSlim Fast Data Loaded\n')

%% Convert spatial units to microns if input units were in pixels.
if Input_space_units == 'pixels'
   for i = 1:size(trajectory,2)
       for j=1:size(trajectory{i},1)
           trajectory{i}{j}(:,2:3) = trajectory{i}{j}(:,2:3)*pixel_size;
       end
   end
end

%From now on, all spatial units are in um.
fprintf('\nPixels converted to microns\n')

%% Filter The Tracks based on their TrackLength.
% Track Length
if TrackMate == 1
   for ff=1:size(trajectory,2)
       for u=1:size(trajectory{1,ff},1)
       Track_length{ff}(u,1) = size(trajectory{1,ff}{u,1},1);
       end
       indices{ff} = find(Track_length{ff}(:,1) >= (minimum_track_length));
       trajectory_filtered{1,ff} = trajectory{1,ff}([indices{ff}]);
   end
end

% Filter tracks from SlimFast based on their tracklength
if SlimFast == 1
   for ff=1:size(trajectory,2)
   indices_length{ff} = find(tracklength{ff} >= (minimum_track_length));
   trajectory_filtered{1,ff} = trajectory{1,ff}([indices_length{ff}]);
   end
end

trajectories = cat(1,trajectory_filtered{:});  %CHOOSE ALL HERE
fprintf('\nTrajectories filtered by length\n')

%% Multiply the frame by the exposure time (in TrackMate, I track them without specifying the frame interval)
if Input_time_units == 'frames'
   for ff=1:size(trajectories,1)
       trajectories{ff,1}(:,1)=trajectories{ff,1}(:,1)*Frame_interval;
   end
   trajectories = trajectories';
else
   trajectories = trajectories';
end
%From now on, trajectories are in um and seconds.
fprintf('\nTrajectories converts to seconds\n')

%% Data Analysis
fprintf('\nCalculating and plotting......\n');


%% ------ Create a separate variable for All the Trajectories -----------
try
ma_AllTracks = ma_AllTracks.addAll(trajectories);
ma_AllTracks = ma_AllTracks.computeMSD;
ma_AllTracks = ma_AllTracks.LogTMSD(TLOGLOG_fitting_points);
ma_AllTracks = ma_AllTracks.TMSD(TMSD_fitting_points);
catch
  fprintf('Error Calculating MSD for All Tracks');
end


%% Identify Butterfly Tracks --------------------------------------
%--------------------------------------------------------------------------------------------
%%1. Detect "butterfly" motion
try
[Butterfly_trajectories,Butterfly_track_segments,Butterfly_Reference_Track_segment, motion_type_segment,reference_of_Tracks_butterfly,Num_BigJumps,BigJump_idx]=identify_Butterfly_tracks_V4(trajectories,jump_threshold,minim_dist,TLOGLOG_fitting_points,R2LIMIT,min_alpha_directed,Offset0,num_points,D0,R0,Conf2JumpThr,Out_percentage,N_sliding,P_min_linear,angleTH,Max_Jump_Confined_Butterfly);
trajectories(reference_of_Tracks_butterfly)=[];

if isempty(Butterfly_trajectories);
    fprintf('There are no Butterfly trajectories');
else
Butterfly_trajectories_segments = horzcat(Butterfly_track_segments{:}); %Reshape cell array to colapse all the cells in the top level.
motion_type_segment = horzcat(motion_type_segment{:}); %Reshape cell array to colapse all the cells in the top level.

%2. Discard those Butterfly tracks with confined segments shorter than N points.
discard_idxxxx = [];
for i=1:length(Butterfly_trajectories);
idxxxx = find(Butterfly_Reference_Track_segment == i);
motion_type_segment_local = motion_type_segment(idxxxx);
idxxxx2 = find(motion_type_segment_local == 1);

if ~isempty(idxxxx2);
   for uu = 1:length(idxxxx2);
       confined_segments_length_local(uu) = length(Butterfly_trajectories_segments{idxxxx(idxxxx2(uu))});
   end
   if ~isempty(confined_segments_length_local);
   if ~isempty(find(confined_segments_length_local < Min_num_points_butt_confined_segment));
       discard_idx(i) = 1;
       discard_idxxxx = [discard_idxxxx idxxxx];
   else
       discard_idx(i) = 0;
   end
   end
end
   clear confined_segments_length_local;
end

Butterfly_trajectories(find(discard_idx==1)) = [];
Butterfly_Reference_Track_segment(discard_idxxxx) = [];
motion_type_segment(discard_idxxxx) = [];
Butterfly_trajectories_segments(discard_idxxxx) = [];


%%3. Add The trajectories to msdanalyzer
ma_butterfly = ma_butterfly.addAll(Butterfly_trajectories);


% FIGURE 1: Plotting of the T-MSD for Butterfly Tracks only
figure()
ma_butterfly.plotMSD;
title([baseName {' '} 'T-MSD Butterfly Tracks']);
filename = strcat(baseName, '_Butterfly_PlotT-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
% xlim([0 xmax]);
% ylim([0 ymax]);


%FIGURE 2: Plotting of the mean T-MSD (black) and standard deviation (Grey) for Butterfly Tracks only
figure()
ma_butterfly.plotMeanMSD(gca, true);
title([baseName {' '} 'TE-MSD Butterfly Tracks']);
filename = strcat(baseName, '_Butterfly_PlotMeanT-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')

% xlim([0 xmax]);
% ylim([0 ymax]);


%FIGURE 3: Visualization of Butterfly Tracks
figure()
ma_butterfly.plotTracks;
ma_butterfly.labelPlotTracks;
title([baseName {' '} 'Butterfly Tracks']);
filename = strcat(baseName, '_Butterfly_PlotTracks');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')

% set(gca,'Ydir','reverse');
%fprintf('Test plotMotion.m')
%plotMotion(ma_butterfly,baseName,OutputFolder)
%fprintf('Works...')
ma_butterfly_segments = ma_butterfly_segments.addAll(Butterfly_trajectories_segments');

%FIGURE 4: Visualization of the confined component of Butterfly Tracks
figure()
plotTracks(ma_butterfly_segments,gca,find(motion_type_segment==1));
ma_butterfly_segments.labelPlotTracks;
title([baseName {''} 'Butterfly Tracks Segmentated (Confined Segments)']);
filename = strcat(baseName, '_Butterfly-confined_PlotTracks');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')

% set(gca,'Ydir','reverse');
%
%
% %FIGURE 5: Visualization of the directed component of Butterfly Tracks
figure()
plotTracks(ma_butterfly_segments,gca,find(motion_type_segment==0));
ma_butterfly_segments.labelPlotTracks;
title([baseName {''} 'Butterfly Tracks Segmentated (Directed Segments)']);
% set(gca,'Ydir','reverse');
filename = strcat(baseName, '_Butterfly-directed_PlotTracks');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')

ma_butterfly_segments_confined = ma_butterfly_segments_confined.addAll(Butterfly_trajectories_segments(find(motion_type_segment==1)));
ma_butterfly_segments_directed = ma_butterfly_segments_directed.addAll(Butterfly_trajectories_segments(find(motion_type_segment==0)));

%Compute TMSD and TEMSD
ma_butterfly = ma_butterfly.computeMSD;
ma_butterfly = ma_butterfly.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly = ma_butterfly.TMSD(TMSD_fitting_points);

ma_butterfly_segments_confined = ma_butterfly_segments_confined.computeMSD;
ma_butterfly_segments_confined = ma_butterfly_segments_confined.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly_segments_confined = ma_butterfly_segments_confined.TMSD(TMSD_fitting_points);

ma_butterfly_segments_directed = ma_butterfly_segments_directed.computeMSD;
ma_butterfly_segments_directed = ma_butterfly_segments_directed.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly_segments_directed = ma_butterfly_segments_directed.TMSD(TMSD_fitting_points);
end
catch
   fprintf('Error with Butterfly Analysis. Likely too few particles');
end

% %% ------------------------------------------------------------------------------
% %Add the trajectories to the msdAnalyzer-----------------------------------
% %------------------------------------------------------------------------------
if isempty(trajectories)==1;
       error('There are no tracks that fulfill the requirements')
else

   ma = ma.addAll(trajectories);
   Num_Tracks = size(trajectories,2);

   for j=1:Num_Tracks;
       Track_Length(j) = size(trajectories{j},1);
       Residence_Times(j)=Track_Length(j)*Frame_interval;
   end

end

% %% -------------------------------------------------------------
% %Calculate the distance travelled by the Tracks-------------------------
% %(This value depends on the Track Length of each Track, but can be informative)--------
% %------------------------------------------------------------------
for ggg = 1:size(trajectories,2);
   points_coord = trajectories{ggg}(:,2:3);
   [max_dist{ggg}, min_dist{ggg}, avg_dist{ggg}] = distance_scatter(points_coord);
end



%% Compute TE-MSD
try
ma = ma.computeMSD;
catch
   fprintf('Error with Computing TE-MSD');
end

%% Preliminary Plotting
% FIGURE 5: Visualization of all the motion tracks
try
figure()
ma_AllTracks.plotTracks;
ma_AllTracks.labelPlotTracks;
set(gca,'Ydir','reverse');
title([baseName {''} 'All Trajectories'])
filename = strcat(baseName, '_AllTracks_PlotTracks');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('All Motion Plot error');
end

% FIGURE 6: Plot of the average T-MSD (black) and standard deviation (grey) for all motion tracks
try
figure()
ma.plotMeanMSD(gca, true)
mmsd = ma.getMeanMSD;
temps = mmsd(:,1);
xs = mmsd(:,2);
dx_plot = mmsd(:,3) ./ sqrt(mmsd(:,4));
dxs = mmsd(:,3);
title([baseName {''} 'All Trajectories'])
filename = strcat(baseName, '_AllTracks_PlotMeanT-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('T-MSD plot error');
end

% %FIGURE 7: Distribution of tracklength for all motion tracks
try
figure()
hist(Track_Length,100);
title([baseName {''} 'Histogram of All Motion Track Lengths']);
xlabel('Track length (frames)');
ylabel('Frequency');
xlim([0 100]);
set(gca,'FontSize',20,'FontWeight','bold');
filename = strcat(baseName, '_AllTracks_TrackLength');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('All motion distribution');
end

%% Compute T-MSD
%This other approach fits every single MSD curve and then the histogram
%of D is very informative
try
ma = ma.TMSD(TMSD_fitting_points);
good_enough_fit_Ds = find(ma.lfit.r2fit >= R2LIMIT);
Dmean = mean( ma.lfit.a(good_enough_fit_Ds) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit_Ds) ) / 2 / ma.n_dim;
fprintf('**Estimation of the diffusion coefficient from linear fit of the MSD curves (Fitting every MSD curve)**:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
   Dmean, Dstd, length(good_enough_fit_Ds));
Ds = ma.lfit.a(good_enough_fit_Ds)/ 2 / ma.n_dim;
catch
   fprintf('Compute T-MSD Error');
end
%
% %% FIGURE 9: Distribution of Diffusion Coefficients for all motion tracks
try
figure()
%Take out negative values from the Diffusion Coefficients list
idx2=find(Ds > 0);
histogram(log10(Ds(idx2)),50);
%xlim([0 0.01]);
title([baseName {''} 'Histogram of diffusion coefficients']);
xlabel('Log10(Diffusion Coefficient)')
ylabel('Frequency');
set(gca,'FontSize',20,'FontWeight','bold');
filename = strcat(baseName, '_AllTracks_DiffCoef');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename))


catch
   fprintf('Distribution of Diff Coeff Error');
end


%% Motion type analysis from T-MSD curves (fitting each Track MSD)
try
ma = ma.LogTMSD(TLOGLOG_fitting_points);

r2fits = ma.loglogfit.r2fit;
alphas = ma.loglogfit.alpha;

% Remove bad fits
bad_fits = r2fits < R2LIMIT;
good_enough_fit_alpha = find(r2fits >= R2LIMIT);
fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
alphas_filtered = alphas(good_enough_fit_alpha);
%Remove NaN from alphas
catch
   fprintf('LogTMSD Error');
end

%% Separate and plot Tracks based on their alpha value from LogLog fit of the T-MSD---------------
%-------------------------------------------------------------------------------------
try
temporal_tray = trajectories(good_enough_fit_alpha);
conf_filtered = find(alphas_filtered <= max_alpha_conf);
brownian_filtered = find(alphas_filtered>max_alpha_conf & alphas_filtered<min_alpha_directed);
directed_filtered = find(alphas_filtered >= min_alpha_directed);
catch
   fprintf('Error with LogLog Fit of MSD');
end
%---------------------------------------------------------------------------------
%------------------------------------------------------------------------------
try
Conf_trajectories_filtered = temporal_tray(conf_filtered);
catch
   fprintf('Calculate # of Tracks per motion type');
end
try
Directed_trajectories_filtered = temporal_tray(directed_filtered);
catch
   fprintf('Calculate # of Tracks per motion type');
end
try
Brownian_trajectories_filtered = temporal_tray(brownian_filtered);
catch
   fprintf('Error with filtered');
end
try
ma_confined = ma_confined.addAll(Conf_trajectories_filtered);
catch
   fprintf('Calculate # of Tracks per motion type');
end
try
ma_directed = ma_directed.addAll(Directed_trajectories_filtered);
catch
   fprintf('Calculate # of Tracks per motion type');
end
try
ma_brownian = ma_brownian.addAll(Brownian_trajectories_filtered);
catch
   fprintf('');
end

%Calculate the number of Tracks per type of motion
try
Percentage_confined = size(ma_confined.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
catch
   fprintf('Calculate # of Tracks per motion type');
end
try
Percentage_brownian = size(ma_brownian.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
catch
   fprintf('Calculate # of Tracks per motion type');
end
try
Percentage_directed = size(ma_directed.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
catch
   fprintf('Calculate # of Tracks per motion type');
end
try
Percentage_butterfly = size(ma_butterfly.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
catch
   fprintf('Calculate # of Tracks per motion type');
end

%Compute MSD for all the motion types in separate variables
try
if isempty(ma_confined.tracks);
else
ma_confined = ma_confined.computeMSD;
ma_confined = ma_confined.LogTMSD(TLOGLOG_fitting_points);
ma_confined = ma_confined.TMSD(TMSD_fitting_points);
end
catch
   fprintf('');
end
try
if isempty(ma_directed.tracks);
else
ma_directed = ma_directed.computeMSD;
ma_directed = ma_directed.LogTMSD(TLOGLOG_fitting_points);
ma_directed = ma_directed.TMSD(TMSD_fitting_points);
end
catch
   fprintf('');
end
try
if isempty(ma_brownian.tracks);
else
ma_brownian= ma_brownian.computeMSD;
ma_brownian = ma_brownian.LogTMSD(TLOGLOG_fitting_points);
ma_brownian = ma_brownian.TMSD(TMSD_fitting_points);
end
catch
   fprintf('');
end

%% Separate the Confined Tracks based on their Mean Jump to eliminate "false confined" tracks.
try
for i=1:length(Conf_trajectories_filtered);
   tracktemp = Conf_trajectories_filtered{i};
   jd_temp_confined_mean(i) = mean(sqrt((tracktemp(2:end,2) - tracktemp(1:end-1,2)).^2 + (tracktemp(2:end,3) - tracktemp(1:end-1,3)).^2));
end
idx6 = find(jd_temp_confined_mean >= Jump_confined_threshold);
idx7 = find(jd_temp_confined_mean < Jump_confined_threshold);

trajectories_confined_High_D = Conf_trajectories_filtered(idx6);
ma_confined_High_D = ma_confined_High_D.addAll(trajectories_confined_High_D);

trajectories_confined_Low_D = Conf_trajectories_filtered(idx7);
ma_confined_Low_D = ma_confined_Low_D.addAll(trajectories_confined_Low_D);


if isempty(ma_confined_High_D.tracks);
else
ma_confined_High_D = ma_confined_High_D.computeMSD;
ma_confined_High_D = ma_confined_High_D.LogTMSD(TLOGLOG_fitting_points);
ma_confined_High_D = ma_confined_High_D.TMSD(TMSD_fitting_points);

%FIGURE 10: Visualize "false confined" motions
% figure()
% ma_confined_High_D.plotTracks;
% ma_confined_High_D.labelPlotTracks;
% title('Confined Trajectories with Average Jump higher than Jump_threshold','Interpreter', 'none');
% set(gca,'Ydir','reverse');


end
catch
   fprintf('');
end
try
if isempty(ma_confined_Low_D.tracks);
else
ma_confined_Low_D = ma_confined_Low_D.computeMSD;
ma_confined_Low_D = ma_confined_Low_D.LogTMSD(TLOGLOG_fitting_points);
ma_confined_Low_D = ma_confined_Low_D.TMSD(TMSD_fitting_points);

%FIGURE 11: Visualize "true confined" tracks
% figure()
% ma_confined_Low_D.plotTracks;
% ma_confined_Low_D.labelPlotTracks;
% title('Confined Trajectories with Average Jump lower than Jump_threshold','Interpreter', 'none');
% set(gca,'Ydir','reverse');
end
catch
   fprintf('');
end
try
if isempty(ma_butterfly_segments.tracks);
else
ma_butterfly_segments = ma_butterfly_segments.computeMSD;
ma_butterfly_segments = ma_butterfly_segments.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly_segments = ma_butterfly_segments.TMSD(TMSD_fitting_points);
end
catch
   fprintf('');
end



%% Plotting ---------------------------------------------------
%--------------------------------------------------------------------------------------
%--------------------------------------------
%-------------------------------------------------------------------------------------
%FIGURE 12: Plot T-MSD of unfiltered confined tracks
try
figure()
ma_confined.plotMSD;
% xlim([0 xmax]);
% ylim([0 ymax]);
title([baseName {''} 'T-MSD Confined Tracks']);
filename = strcat(baseName, '_Confined_DiffCoef');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('');
end

% %FIGURE 12: Plot average T-MSD (black) and standard deviation (grey) of unfiltered confined tracks
try
figure()
ma_confined.plotMeanMSD(gca, true);
% xlim([0 xmax]);
% ylim([0 ymax]);
title([baseName {''} 'TE-MSD Confined Tracks']);
filename = strcat(baseName, '_confined_PlotTE-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('');
end

%FIGURE 13: Visualize unfiltered confined tracks
try
figure()
ma_confined.plotTracks;
ma_confined.labelPlotTracks;
title('Confined Tracks');
set(gca,'Ydir','reverse');
catch
   fprintf('');
end
%FIGURE 14: Plot T-MSD for directed motion tracks
try
figure()
ma_directed.plotMSD;
% xlim([0 xmax]);
% ylim([0 ymax]);
title([baseName {''} 'T-MSD Directed Tracks']);
filename = strcat(baseName, '_directed_PlotT-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('');
end

%FIGURE 15: Plot average T-MSD (black) and standard deviation (grey) of
% directed motion tracks
try
figure()
ma_directed.plotMeanMSD(gca, true);
% xlim([0 xmax]);
% ylim([0 ymax]);
title([baseName {''} 'TE-MSD Directed Tracks']);
filename = strcat(baseName, '_directed_PlotTE-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('');
end

%FIGURE 16: Visualize directed motion tracks
try
figure()
ma_directed.plotTracks;
ma_directed.labelPlotTracks;
title([baseName {''} 'Directed Tracks']);
set(gca,'Ydir','reverse');
filename = strcat(baseName, '_directed_PlotTracks');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('');
end
%FIGURE 17: Plot T-MSD for pure Brownian motion tracks
try
figure()
ma_brownian.plotMSD;
title([baseName {''} 'T-MSD Brownian Tracks']);
filename = strcat(baseName, '_Brownian_PlotT-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
% xlim([0 xmax]);
% ylim([0 ymax]);
catch
   fprintf('');
end
%FIGURE 18: Plot average T-MSD (black) and standard deviation (grey) for pure Brownian motion tracks
try
figure()
ma_brownian.plotMeanMSD(gca, true);
title([baseName {''} 'TE-MSD Brownian Tracks']);
filename = strcat(baseName, '_Brownian_PlotTE-MSD');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
% xlim([0 xmax]);
% ylim([0 ymax]);
catch
   fprintf('');
end
%FIGURE 19: Visualize pure Brownian tracks
try
figure()
ma_brownian.plotTracks;
ma_brownian.labelPlotTracks;
title([baseName {''} 'Brownian Tracks']);
set(gca,'Ydir','reverse');
filename = strcat(baseName, '_Brownian_PlotT');
saveas(gcf, fullfile(OutputFolder,filename))
saveas(gcf, fullfile(OutputFolder,filename), 'tif')
catch
   fprintf('');
end


%%   save the results from the MSD Analysis
fprintf('\nSaving MSD Results...\n');
clear ma

% Save Input Variables
FileName = fullfile(OutputFolder,strcat(baseName,'_InputVariables','.csv'));
writetable(struct2table(inputVar), FileName)

% ALL TRACKS
try
ma = ma_AllTracks;
fileName = strcat(baseName,'_msd_results_','AllTracks','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving All Tracks');
end
% BROWNIAN
try
ma = ma_brownian;
fileName = strcat(baseName,'_msd_results_','Brownian','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving Brownian Tracks');
end
% CONFINED
try
ma = ma_confined;
fileName = strcat(baseName,'_msd_results_','Confined','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving Confined Tracks');
end
% DIRECTED
try
ma = ma_directed;
fileName = strcat(baseName,'_msd_results_','Directed','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving Directed Tracks');
end

% BUTTERFLY
try
ma = ma_butterfly;
fileName = strcat(baseName,'_msd_results_','Butterfly','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving Butterfly Tracks');
end
% CONFINED - HIGH DIFF
try
ma = ma_confined_High_D;
fileName = strcat(baseName,'_msd_results_','Confined_HighDif','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving Confined - High Dif');
end
% CONFINED - LOW DIFF
try
ma = ma_confined_Low_D;
fileName = strcat(baseName,'_msd_results_','Confined_LowDif','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving Confined - Low Dif');
end
% BUTTERFLY - CONFINED SEG
try
ma = ma_butterfly_segments_confined;
fileName = strcat(baseName,'_msd_results_','Butterfly_Confined_Segments','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear ma fileName;
catch
   fprintf('Error Saving Butterfly Confined Segments');
end
% BUTTERFLY - DIRECTED SEG
try
ma = ma_butterfly_segments_directed;
fileName = strcat(baseName,'_msd_results_','Butterfly_Directed_Segments','.mat');
save(fullfile(OutputFolder,fileName),'ma');
clear fileName
catch
   fprintf('Error Saving Butterfly Directed Segments');
end
% Butterfly Segments - Ref ID
try
fileName = strcat(baseName,'_msd_results_','Butterfly_Segments_RefID','.mat');
save(fullfile(OutputFolder,fileName),'ma');
catch
   fprintf('Error Saving Butterfly Ref IDs');
end
end
```
**MSD_LinearFit.m**
```Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ALL TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Saving Information
  % settings
    %directory_name = uigetdir('C:\Users\letinunez\Documents\');% select output folder
    % directory_name ='/Users/letinunez/ACTBmRNAZBP1KOneurons/inst/Analysis_28APR2022'; %MAC
    directory_name ='E:\New Dropbox\Dropbox (EinsteinMed)\ACTBmRNAZBP1KOneurons\inst\Analysis_29APR2022'; % PC- HOME
    date = datestr(datetime, 'yyyy-mm-dd_HH-MM');

  % Folder and File Names
    folderName = fullfile(directory_name,strcat('Analyzed_Results_LinearFit_',date));
  % File name for all 'ma' files to be stored (e.g. ma_AllTracks,
  % ma_directed)
   % filebaseName =  strcat('msd_results_',num2str(Frame_interval),'ms_',num2str(Jump_confined_threshold),'um-jumpthres_', num2str(R2LIMIT),'-R2Limit');
  % make folder in directory
    mkdir(folderName);
%% Plot MSD Results for WT and KO on same graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1. All Tracks

% set data type set
%motion = 1;
% 1 - AllTracks
% 2 - Confined 
% 3 - Confined - H
% 4 - Confined - L
% 5 - Directed
% 6 - Brownian
% 7 - Butterfly
% 8 - Butterfly - Confined
% 9 - Butterfly - Directed
%Set the folder as a cell array. (If only one
% location is included, then it will be used for all data types)


% Inputs
%------------------------------------------------------------------
%Data---------------------------------------------------
% %------------------------------------------------------------------
%Set the data type description (Maximum 8 types)
type.t1 = 'WT'; % condition 1
type.t2 = 'KO'; % condition 2



% startlocn = uigetdir('C:\Users\letinunez\Documents\'); % save selected files together as a group
% startlocn = '/Users/letinunez/ACTBmRNAZBP1KOneurons/inst/Analysis_24APR2022/MSD_Results';
startlocn = 'E:\New Dropbox\Dropbox (EinsteinMed)\ACTBmRNAZBP1KOneurons\inst\Analysis_29APR2022\MSD_Results'; % PC- HOME
%Save_Results = 0; % 1 - yes save results; 0 - no
%File_name = 'msd_results'; % prefix to msd .csv files


%% 
fprintf('Initializing Input Variables\n');
%------------------------------------------------------------------
%Numerical Inputs---------------------------------------------------
%------------------------------------------------------------------
n_dim = 2; %Dimensionality of the movement (2 for 2D and 3 for 3D).
Frame_interval = 0.10; %Exposure time in seconds

TMSD_fitting_points = 3;   %Minimum 3 point for being able to calculate the Confidence Intervals
TEMSD_fitting_points = 3;  %Minimum 3 point for being able to calculate the Confidence Intervals
TLOGLOG_fitting_points = 20; 

R2LIMIT = 0.7; %R-squared minimum value from the D fits on each T-MSD Curve.



%-----------------------------------------------------------------
%For the Confined Circle Diffusion Model Fitting---------------------------------------------------
%---------------------------------------------------------------------
num_points = 10; %Number of points for fitting the confined diffussion circle model to the TE-MSD.
level = 0.001; %This is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in ?m]. 
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (?m^2/s)
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (?m)

Radius_interval = [30 60]; %Radius interval (in nm) in order to divide tracks into several populations. (log10(60) = 1.78)



%------------------------------------------------------------------------
%Initialize variables-------------------------------------------------
%----------------------------------------------------------------------
% Dcoeffs = {};
datatypes = fieldnames(type);
ntypes = size(datatypes,1);    
for i=1:ntypes;
files{i}= Select1DataGroup(type.(datatypes{i}),'mat',startlocn);
end

%n = length(MSD_Results_total.ma_AllTracks);
SPACE_UNITS = 'µm'; %This is just for visualization purposes
TIME_UNITS = 's'; %This is just for visualization purposes
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
fprintf('Starting MSD Linear Fit\n');

for i=1:ntypes;
MSD_Results{i}.ma = msdanalyzer(2, 'um', 's');  %Initialize MSD analizer
f = waitbar(0, 'Starting');
n = size(files{i}.data,1);
    for c=1:n;
        MSD_Results_prov{i}{c} = load(strcat(files{i}.data{c,2},files{i}.data{c,1}));
        [folder, baseName,ext] = fileparts(files{1,i}.data{c,1});
        fprintf(1, 'Now reading %s\n', baseName);

        ma = MSD_Results_prov{1,i}{1,c}.ma;
        figure
        try
        ma.plotTracks
        ma.labelPlotTracks
        title(baseName);
        filename = strcat(baseName, '_PlotTracks');
        saveas(gcf, fullfile(folderName,filename))
        saveas(gcf, fullfile(folderName,filename), 'tif')

        figure
        ma.plotMSD
        title(baseName);
        filename = strcat(baseName, '_PlotMSD');
        saveas(gcf, fullfile(folderName,filename))
        saveas(gcf, fullfile(folderName,filename), 'tif')

        figure
        ma.plotMeanMSD(gca, true)
        title(baseName);
        % Adds linear regression to MeanMSD plot
        %[fo, gof] = ma.fitMeanMSD( 0.25 );
        %plot(fo)
        filename = strcat(baseName, '_PlotMeanMSD');
        saveas(gcf, fullfile(folderName,filename))
        saveas(gcf, fullfile(folderName,filename), 'tif')
        
        % Linear Fit
        ma = MSD_Results_prov{1,i}{1,c}.ma.TMSD(TMSD_fitting_points);
        good_enough_fit_Ds_2{1,i} = find(ma.lfit.r2fit >= R2LIMIT);
        Dmean_2{1,i} = mean( ma.lfit.a(good_enough_fit_Ds_2{1,i}) ) / 2 / ma.n_dim;
        Dstd_2{1,i}  =  std( ma.lfit.a(good_enough_fit_Ds_2{1,i}) ) / 2 / ma.n_dim;
        fprintf('**Estimation of the diffusion coefficient from linear fit of the MSD curves (Fitting every MSD curve)**:\n')
        fprintf('D = %.3g +/- %.3g (mean +/- std, N = %d)\n', ...
        Dmean_2{1,i}, Dstd_2{1,i}, length(good_enough_fit_Ds_2{1,i}));
        Ds_2{1,i} = ma.lfit.a(good_enough_fit_Ds_2{1,i})/ 2 / ma.n_dim;
        
        filename = strcat(baseName,'_DiffusionCoeff.csv');
        dlmwrite(fullfile(folderName,filename),Ds_2{1,i})
        catch
            fprintf('Error - Too few Tracks');
        end
        % Progress Bar
        waitbar(c/n, f, sprintf('%s Progress: %d %%',type.(datatypes{i}), floor(c/n*100)));
        pause(0.1);
        close all
        clearvars -except folderName MSD_Results files i c SPACE_UNITS TIME_UNITS ma datatypes ntypes TMSD_fitting_points R2LIMIT n_dim n f type
     
    end
    close(f) %closes progress bar
end
fprintf('Finished MSD Linear Fit\n');
%% All Tracks - Histogram of Diffusion Coefficients - WORKING


% figure()
% %Take out negative values from the Diffusion Coefficients list
% idx2_2=find(Ds_2{1,1} > 0);
% histogram(log10(Ds_2{1,1}(idx2_2)),50);
% %xlim([0 0.01]);
% hold on
% idx2_2=find(Ds_2{1,2} > 0);
% histogram(log10(Ds_2{1,2}(idx2_2)),50);
% hold off
% trackname = 'All Tracks';
% title([trackname {''} 'Histogram of diffusion coefficients']);
% xlabel('Log10(Diffusion Coefficient)')
% ylabel('Frequency');
% legend('WT','KO');
% name = 'combined';
% filename = strcat(name, '_DcoefHisto');
% saveas(gcf, fullfile(folderName,filename))
% saveas(gcf, fullfile(folderName,filename), 'tif')
% %set(gca,'FontSize',20,'FontWeight','bold');
%%
%filename = strcat(name,'_DiffusionCoeff_AllTracks.csv');
%dlmwrite(fullfile(folderName,filename),log10(Ds_2{1,2}(idx2_2)))

%% Ideallyy export data to csv so that there can be one file with the diffusion coefficients
% check if file exists
```
**MSD_ParabolicFit.m**
```Matlab
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ALL TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Saving Information
  % settings
  %  directory_name = uigetdir('C:\Users\letinunez\Documents\');% select output folder
  %  directory_name ='/Users/letinunez/ACTBmRNAZBP1KOneurons/inst/Analysis_24APR2022'; %MAC
    directory_name ='E:\New Dropbox\Dropbox (EinsteinMed)\ACTBmRNAZBP1KOneurons\inst\Analysis_29APR2022'; % PC- HOME

    date = datestr(datetime, 'yyyy-mm-dd_HH-MM');

  % Folder and File Names
    folderName = fullfile(directory_name,strcat('Analyzed_Results_ParabolicFit_',date));
  % File name for all 'ma' files to be stored (e.g. ma_AllTracks,
  % ma_directed)
   % filebaseName =  strcat('msd_results_',num2str(Frame_interval),'ms_',num2str(Jump_confined_threshold),'um-jumpthres_', num2str(R2LIMIT),'-R2Limit');
  % make folder in directory
    mkdir(folderName);
%% Plot MSD Results for WT and KO on same graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1. All Tracks

% set data type set
%motion = 1;
% 1 - AllTracks
% 2 - Confined 
% 3 - Confined - H
% 4 - Confined - L
% 5 - Directed
% 6 - Brownian
% 7 - Butterfly
% 8 - Butterfly - Confined
% 9 - Butterfly - Directed
%Set the folder as a cell array. (If only one
% location is included, then it will be used for all data types)


% Inputs
%------------------------------------------------------------------
%Data---------------------------------------------------
% %------------------------------------------------------------------
%Set the data type description (Maximum 8 types)
type.t1 = 'WT'; % condition 1
type.t2 = 'KO'; % condition 2



%startlocn = uigetdir('C:\Users\letinunez\Documents\'); % save selected files together as a group
startlocn = '/Users/letinunez/ACTBmRNAZBP1KOneurons/inst/Analysis_24APR2022/MSD_Results';
%Save_Results = 0; % 1 - yes save results; 0 - no
%File_name = 'msd_results'; % prefix to msd .csv files


%% 

%------------------------------------------------------------------
%Numerical Inputs---------------------------------------------------
%------------------------------------------------------------------
n_dim = 2; %Dimensionality of the movement (2 for 2D and 3 for 3D).
Frame_interval = 0.10; %Exposure time in seconds

TMSD_fitting_points = 3;   %Minimum 3 point for being able to calculate the Confidence Intervals
TEMSD_fitting_points = 3;  %Minimum 3 point for being able to calculate the Confidence Intervals
TLOGLOG_fitting_points = 20; 

R2LIMIT = 0.7; %R-squared minimum value from the D fits on each T-MSD Curve.



%-----------------------------------------------------------------
%For the Confined Circle Diffusion Model Fitting---------------------------------------------------
%---------------------------------------------------------------------
num_points = 10; %Number of points for fitting the confined diffussion circle model to the TE-MSD.
level = 0.001; %This is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in ?m]. 
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (?m^2/s)
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (?m)

Radius_interval = [30 60]; %Radius interval (in nm) in order to divide tracks into several populations. (log10(60) = 1.78)



%------------------------------------------------------------------------
%Initialize variables-------------------------------------------------
%----------------------------------------------------------------------
% Dcoeffs = {};
datatypes = fieldnames(type);
ntypes = size(datatypes,1);    
for i=1:ntypes;
files{i}= Select1DataGroup(type.(datatypes{i}),'mat',startlocn);
end

%n = length(MSD_Results_total.ma_AllTracks);
SPACE_UNITS = 'µm'; %This is just for visualization purposes
TIME_UNITS = 's'; %This is just for visualization purposes
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
fprintf('Starting MSD Parabolic Fit\n');

for i=1:ntypes;
MSD_Results{i}.ma = msdanalyzer(2, 'um', 's');  %Initialize MSD analizer
f = waitbar(0, 'Starting');
n = size(files{i}.data,1);
    for c=1:n;
        MSD_Results_prov{i}{c} = load(strcat(files{i}.data{c,2},files{i}.data{c,1}));
        [folder, baseName,ext] = fileparts(files{1,i}.data{c,1});
        ma = MSD_Results_prov{1,i}{1,c}.ma;
        fprintf(strcat('\nParabolic Fit for: ',baseName,'\n'));
%         figure
        try
%         ma.plotTracks
%         ma.labelPlotTracks
%         title(baseName);
%         filename = strcat(baseName, '_PlotTracks');
%         saveas(gcf, fullfile(folderName,filename))
%         saveas(gcf, fullfile(folderName,filename), 'tif')
% 
%         figure
%         ma.plotMSD
%         title(baseName);
%         filename = strcat(baseName, '_PlotMSD');
%         saveas(gcf, fullfile(folderName,filename))
%         saveas(gcf, fullfile(folderName,filename), 'tif')
% 
%         figure
%         ma.plotMeanMSD(gca, true)
%         title(baseName);
%         filename = strcat(baseName, '_PlotMeanMSD');
%         saveas(gcf, fullfile(folderName,filename))
%         saveas(gcf, fullfile(folderName,filename), 'tif')
        
        %% Add Parabolic Fitting to Directed Motion Tracks
        A = rmmissing(ma.getMeanMSD); % get Mean MSD and remove missing values (e.g. NaN)
        SPACE_UNITS = 'µm'; %This is just for visualization purposes
        TIME_UNITS = 's'; %This is just for visualization purposes
    
        t = A(:, 1); % delay vector
        msd = A(:,2); % msd
        std_msd = A(:,3); % we will use inverse of the std as weights for the fit
        std_msd(1) = std_msd(2); % avoid infinity weight

        figure
        ma.plotMeanMSD(gca, true)

        ft = fittype('a*x + c*x^2'); % fit to parabola
        [fo, gof] = fit(t, msd, ft, 'Weights', 1./std_msd, 'StartPoint', [0 0],'Exclude', t>0.5);


        hold on
        plot(fo)
        legend off
        ma.labelPlotMSD
        title(baseName);
        filename = strcat(baseName, '_PlotMeanMSD_parabolicfit');
        saveas(gcf, fullfile(folderName,filename))
        saveas(gcf, fullfile(folderName,filename), 'tif')

        Dfit = fo.a / 4;
        Vfit = sqrt(fo.c);

        ci = confint(fo);
        Dci = ci(:,1) / 4;
        Vci = sqrt(ci(:,2));

        fprintf('Parabolic fit of the average MSD curve with 95% confidence interval:\n')

        % Values from MSD Analyzer tutorial
        % Diffusion coefficient. Will set the amplitude of the random displacement
        D  = 1e-3; % µm^2/s
        % Mean velocity
        vm = 0.05; % µm/s

        fprintf('\nD = %.3g [ %.3g - %.3g ] %s, random displacement is %.3g %s\n', ...
        Dfit, Dci(1), Dci(2), [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);

        fprintf('V = %.3g [ %.3g - %.3g ] %s, random velocity is %.3g %s\n', ...
        Vfit, Vci(1), Vci(2), [SPACE_UNITS '/' TIME_UNITS], vm, [SPACE_UNITS '/' TIME_UNITS]);
 
       
        
        filename = strcat(baseName,'_DiffusionCoeff.csv');
        dlmwrite(fullfile(folderName,filename),Dfit)
        catch
            fprintf('Error - Too few Tracks');
        end
        waitbar(c/n, f, sprintf('%s Progress: %d %%',type.(datatypes{i}), floor(c/n*100)));
        pause(0.1);
        close all
        clearvars -except folderName MSD_Results files i c SPACE_UNITS TIME_UNITS ma datatypes ntypes TMSD_fitting_points R2LIMIT n_dim Dfit n f type
    end
    close(f) %closes progress bar
end
fprintf('Finished MSD Parabolic Fit\n');

%% All Tracks - Histogram of Diffusion Coefficients - WORKING


% figure()
% %Take out negative values from the Diffusion Coefficients list
% idx2_2=find(Ds_2{1,1} > 0);
% histogram(log10(Ds_2{1,1}(idx2_2)),50);
% %xlim([0 0.01]);
% hold on
% idx2_2=find(Ds_2{1,2} > 0);
% histogram(log10(Ds_2{1,2}(idx2_2)),50);
% hold off
% trackname = 'All Tracks';
% title([trackname {''} 'Histogram of diffusion coefficients']);
% xlabel('Log10(Diffusion Coefficient)')
% ylabel('Frequency');
% legend('WT','KO');
% name = 'combined';
% filename = strcat(name, '_DcoefHisto');
% saveas(gcf, fullfile(folderName,filename))
% saveas(gcf, fullfile(folderName,filename), 'tif')
% %set(gca,'FontSize',20,'FontWeight','bold');
%%
%filename = strcat(name,'_DiffusionCoeff_AllTracks.csv');
%dlmwrite(fullfile(folderName,filename),log10(Ds_2{1,2}(idx2_2)))

%% Ideallyy export data to csv so that there can be one file with the diffusion coefficients
% check if file exists
```
