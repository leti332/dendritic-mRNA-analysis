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

