% controller
function z = MC_segm(folder, firstFrame, roi, ~)

% Set up parameters to iterate though frames
filelist = dir(folder);
nFiles = size(filelist,1)-2;
prevFrame = strcat([folder, firstFrame]);
prevIm = imread(prevFrame);

% Instantiate the object image

imageObj = ImageClass(folder, 1, prevIm, roi);  

% Begin at 4 because we have a subfolder GT.
for fileNumber = 4:1:nFiles
    tic
    display(['processing file ' int2str(fileNumber)-3]);
 
    frameName = filelist(fileNumber).name;
    currIm = imread(strcat([folder, frameName]));
    
    imageObj.fileNumber = fileNumber;
    imageObj.image = currIm;
    
    getSegmentationC1(imageObj);
    
    %imwrite(imageObj.outMask, sprintf('mask_%d.jpg',fileNumber));
    prevIm = currIm;
    toc
    %TODO:write file of mask and segmented object
    %TODO:compute measures for each file
    
end
%DONE:in the loop call segmentation,
%actually we don't need to roi
%TODO:write a csv file with all measures
end
