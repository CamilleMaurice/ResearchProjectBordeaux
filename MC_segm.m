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
for fileNumber = 4:1:nFiles
    frameName = filelist(fileNumber).name;
    currFrame = strcat([folder, frameName]);
    
    labels = getSegmentation(currFrame, prevFrame, roi);
    %TO DO: create a function from labels to roi, to update the roi.
    prevFrame = currFrame;

end
%DONE:in the loop call segmentation,
%TODO:get roi
%TODO:write file of segmented object
%TODO:compute measures for each file
%TODO:write a csv file with all measures




end