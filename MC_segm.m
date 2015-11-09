% controller
function z = MC_segm(folder,first_frame,roi,method)
disp('beginning initialization');

folder
first_frame

% TODO: iterate through all frames
dname = first_frame;
filelist = dir([fileparts(dname) filesep '*.bmp'])
% fileNames = {filelist.name}';
% num_frames = (numel(filelist));
% 
% for i=1:num_frames
%     I = imread(fullfile(data.pathName, fileNames{i})); %to show the first image in the selected folder
%     imshow(I, []);
%     pause(1)
% end

%TODO:in the loop call segmentation,
%TODO:get roi
%TODO:write file of segmented object
%TODO:compute measures for each file
%TODO:write a csv file with all measures




end