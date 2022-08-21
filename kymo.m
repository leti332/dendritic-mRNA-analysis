function kymo(InputPath, filename, OutputPath)
%%
% fprintf(['InputPath: ' InputPath '\n'])
% fprintf(['File name: ' filename '\n'])
% Concatenate InputPath and filename to give imread full path of file
%pathfileName = [InputPath,filename]; % MAC
%pathfileName = [InputPath,'\\', filename];% PC - HOME
pathfileName = fullfile(InputPath, filename); 
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
imwrite(Kymo, fullfile(OutputPath,strcat(baseFileName, '_kymo.tif')), 'tif', 'WriteMode', 'overwrite', 'Compression', 'none'); % PC
% Save kymograph with t and x units (_kymo_units.tif)
%saveas(gcf, [OutputPath, baseFileName, '_kymo_units.tif'])%MAC
saveas(gcf, fullfile(OutputPath,strcat(baseFileName, '_kymo_units.tif'))) %PC
% Save matlab figure (_kymo.fig)- NOT WORKING
%saveas(gcf, [OutputPath, baseFileName, '_kymo.fig'], '.fig')
end