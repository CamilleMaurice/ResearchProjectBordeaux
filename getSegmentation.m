function binarySegmentationMask = getSegmentation(filename,roi)
n_bins = 16;
max_displacement = 2;
window_omega = 6;
std = 0.5;
nLabels = 50;

gaussian = fspecial('gaussian', window_omega, std);

% addpath(genpath('maxflow-v3.0'));

current_frame = imread( filename );
object_region = imcrop( current_frame, roi );
bkg_region = getBkgRegion( roi, current_frame );

[height,width,~] = size( current_frame );
nPixels = height*width;

%LABELS DEFINITION
index = createIndex();
labelCost = createLabelCost(index);

class = zeros (1,nPixels);
labels = zeros(height, width, 3);

%O for bkg 1 for object
segmentation_labels = labels(:,:,1);

displacement_labels = labels(:,:,2:3);

%DATA TERM
%
%TODO:APPEARANCE SIMILARITY AS A FUNCTION OF the displacement
app_simil = zeros ( height, width, (2*max_displacement+1)^2 );
cpt = 1;
for y = -max_displacement:max_displacement
    for x = -max_displacement:max_displacement
        app_simil(:,:,cpt) = getApperanceSimilarity( patch, displaced_patch );
        cpt = cpt +1;
    end
    
end
%Need a function to create window around each pixel of the image.


%APPEARANCE MODEL
%just check format with gc mex
histoBkg = histo3D( bkg_region, n_bins );
histoObj = histo3D( reshape( [object_region], [], 3 ), n_bins);

%Compute Probabilities - normalize histograms
probsObj = histoObj/size(reshape( [object_region], [], 3 ),1);
probsBkg = histoBkg/size(bkg_region,1);

U_fg = zeros(height, width, 1);
U_bg = zeros(height, width, 1);
%OPTIMIZE THIS
% Rcolor = current_frame(:,:,1);cpde appearan
% Gcolor = current_frame(:,:,2);
% Bcolor = current_frame(:,:,3);
for j = 1 : height
    for i = 1 : width
        Rcolor = floor(current_frame(j,i,1)/n_bins)+1;
        Gcolor = floor(current_frame(j,i,2)/n_bins)+1;
        Bcolor = floor(current_frame(j,i,3)/n_bins)+1;
        U_fg(j,i) = -log (probsObj(Rcolor,Gcolor,Bcolor));  
        U_bg(j,i) = -log (probsBkg(Rcolor,Gcolor,Bcolor));  
    end
end
%%APPEARANCE MODEL DONE



%TODO:SMOOTHNESS TERM spatial

%TODO:MOTION COHERENCE
%TODO:ATTRIBUTE COHERENCE

%TODO: SAME WITH TEMPORAL NEIGHBORS
%TODO:MOITON COHERENCE
%TODO: ATTRIBUTE COHERENCE



reshape( [current_frame], [], 1 );

%TODO CALL GC MEX
%gfet the labels
%TODO: construct the binary segmentation mask from the output labels

end

%returns the bckg pixels colors in a 1D vector
function concatenatedImage = getBkgRegion(roi, image)
%This function has been unary tested
%Size of the image = size of the backgorund + size of the object
size_im = size(image);
RoiULx = roi(1);
RoiULy = roi(2);
RoiWidth = roi(3);
RoiHeight = roi(4);

BkgL = imcrop(image,[1, 1, RoiULx-1, size_im(1)]);
BkgU = imcrop(image,[RoiULx, 1, RoiWidth-1, RoiULy-1]);
BkgD = imcrop(image,[RoiULx, RoiULy + RoiHeight + 1, RoiWidth-1, size_im(1)- (RoiULy + RoiHeight +1)]);
BkgR = imcrop(image,[RoiULx+RoiWidth+1, 1, size_im(2)-(RoiULx + RoiWidth +1), size_im(1)]);

%from matrix to vector;
BkgL = reshape([BkgL], [], 3);
BkgU = reshape([BkgU], [], 3);
BkgD = reshape([BkgD], [], 3);
BkgR = reshape([BkgR], [], 3);

concatenatedImage = [BkgL; BkgU; BkgD ; BkgR] ;
end

function index = createIndex()
    index = zeros(50,4);
    for i = 1:1:50
        index(i,1) = i;
        %0 :: background, 1 :: object
        if ( i < 26 )
            index(i,2) = 0;
        else
            index(i,2) = 1;
        end
    end
    dx = -2;
    for i = 1:5:50
        if ( i == 26 )
            dx = -2;
        end
        cpt = 0;
       for dy = -2:1:2
        index(i+cpt, 4) = dy;
        cpt =cpt +1;
       end
       
       index(i:i+5, 3) = dx; 
       dx = dx +1;
    end
 index = index(1:end -1, :);
end

% For initialization all dx, dy are set to 0, we choose the label
% corresponding to bg/fg with displacement 00;
function class = createClass (image, roi, index)
    [h, w, ~] = size(image);
    class = zeros (h, w);
    labelFg00 = getLabel(index, 1, 0, 0);
    labelBg00 = getLabel(index, 0, 0, 0);
    
    for y = 1:h
        for x = 1:w
            if (isInRoi(y, x, roi))
                class(y,x) = labelFg00;
            else
                class(y,x) = labelBg00;
            end
        end
    end
    
    %/!\ the vector cuts the image column wise and not row wize 
    %Further check is needed to be sure it is according to GCMEX
    %/!\
    class = reshape (class, [], 1)';
end

%TODO : Make further tests (especially about < ou <=)
function  [bool]  = isInRoi (y, x, roi)
    RoiULx = roi(1);
    RoiULy = roi(2);
    RoiWidth = roi(3);
    RoiHeight = roi(4);
    
    if ( x >= RoiULx && x < RoiULx + RoiWidth )
        if (y >= RoiULy && y < RoiULy + RoiHeight)    
            bool = true(1);
        else
            bool = false(1);
        end
    else
        bool = false (1);
    end

end

function label = getLabel (index, seg, dx, dy)
%seg must be 0 or 1, dx and dy  btw -2:2 
%add check here
    i4 = find(index(:, 4) == dy);
    i3 = find(index(:, 3) == dx);
    i2 = find(index(:, 2) == seg);
    label = intersect(i4, i3);
    label = intersect(label, i2);    
end


%Create a label cost matrix according to the format wanted by GCMEX.
function labelCost = createLabelCost (index)
    [nLabels, ~] = size (index); 
    labelCost = zeros (nLabels, nLabels);
    for lp = 1:nLabels
        for lq = 1:nLabels
            labelCost(lp,lq) = getDistanceBtwLabels(lp, lq, index);
        end
    end
    
end

%Distance is defined as the number of different informations between 2
%labels. Max distance = 3; Min distance = 0 (if lp == lq)
function distance = getDistanceBtwLabels ( lp, lq, index )
    lp_info = index(lp,2:end);
    lq_info = index(lq,2:end);
    distance = sum(abs(lp_info - lq_info)); 
end


function result = SSDG (curr_patch, displaced_patch, gaussian)


end