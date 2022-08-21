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
