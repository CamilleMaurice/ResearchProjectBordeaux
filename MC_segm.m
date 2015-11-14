% controller
function z = MC_segm(folder,first_frame,roi,method)
disp('beginning initialization');

% DONE: iterate through all frames
dname = first_frame;
filelist = dir(folder);
%For more genericity we should find a way to get rid of this -3, at least
%back to -2 (detect the folder).
nFiles = size(filelist,1)-2;
filelist.name;
prevFrame = strcat([folder, first_frame]);
prev_Im=imread(prevFrame);
[h,w,~]=size(prev_Im)
for fileNumber = 4:1:nFiles%WHY: the beginning at 4?
    display(['processing file ' int2str(fileNumber)]);
    frameName = filelist(fileNumber).name;
    currFrame = strcat([folder, frameName]);
    curr_Im=imread(currFrame);
    
    [labels, binaryMask ] = getSegmentation(curr_Im, prev_Im, roi);
    
    prev_Im = curr_Im;
 
    %TODO:write file of mask and segmented object
    %TODO:compute measures for each file

    
end
%DONE:in the loop call segmentation,
%actually we don't need to roi
%TODO:write a csv file with all measures
end
