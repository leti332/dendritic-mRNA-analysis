%% Creates Average Z-projection of movie
function OutputImage=avg_proj(ZstackNum,output_PathName,PathFileName,option)  % StartImageNum >10,  EndImageNum<100
%% Define Path and File Name
[pathName, fileName] = fileparts(PathFileName); % get path and filename (without extension) separately
%PathName=[pathName,'\']; %MAC
%PathName = pathName, '\\']; % PC
%FileName=[fileName,'.tif'];
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
    imwrite(proj, fullfile(output_PathName,strcat('AVG_',fileName, '.tif')), 'tif', 'WriteMode', 'append', 'Compression', 'none');%PC
end