%% Clear All Variables before Beginning
clear all % clear all variables from workspace
close all % close all figures
clc % clear command window/console
%% Select Directory (Folder)≈
% directs GUI to this folder first
rootFolder =  '/Users/letinunez/Desktop/'; % MAC
%'C:\\Users\\letinunez\\Desktop\\'; % PC
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
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Runs Kymo function to generate Kymograph
  kymo(myDir, baseFileName, OutputPath)
  % Progress Bar
  waitbar(k/n, f, sprintf('Progress: %d %%', floor(k/n*100)));
  pause(0.1);
end
close(f) %closes progress bar
fprintf('Finished \n')