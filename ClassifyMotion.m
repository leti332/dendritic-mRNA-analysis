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
