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
[h,w]=size(prev_Im);
for fileNumber = 4:1:nFiles%WHY: the beginning at 4?
    frameName = filelist(fileNumber).name;
    currFrame = strcat([folder, frameName]);
    curr_Im=imread(currFrame);
    
    %labels = getSegmentation(curr_Im, prev_Im, roi);
    %for testing
    labels = zeros(h,w);
    labels(2:4,2:10)=1;
    
    prev_Im = curr_Im;
    
    mask=fromLabelstoMask(labels, w, h);
    imshow(mask)
    roi =  fromMaskToRoi ( labels, w );
    %TODO:write file of mask and segmented object
    %TODO:compute measures for each file

    
end
%DONE:in the loop call segmentation,
%TODO:get roi
%TODO:write a csv file with all measures
end

%TO DO: create a function from labels to roi, to update the roi.
function roi = fromMaskToRoi ( mask, width )
CC=bwconncomp(mask);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
[col row]=ind2sub(width,CC.PixelIdxList{idx});
bb=regionprops(CC,'BoundingBox');
roi=[bb.BoundingBox(1:2)+0.5 bb.BoundingBox(3:4)];
%TO DO : get the UL coordinate of the obj.
%Find the width of the obj
%Find the height

end

function mask = fromLabelstoMask(labels,width,height)
mask=zeros(height,width);
med=max(max(labels))/2;
mask(labels>med)=1;%check the strict superiority

end