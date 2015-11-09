function binarySegmentationMask = getSegmentation(filename,roi)
n_bins = 16;
max_displacement = 2;
window_omega = 6;
std = 0.5;

gaussian = fspecial('gaussian', window_omega, std);

% addpath(genpath('maxflow-v3.0'));

current_frame = imread( filename );
object_region = imcrop( current_frame, roi );
bkg_region = getBkgRegion( roi, current_frame );

[height,width,~] = size( current_frame );
number_of_pixels = height*width;

%TODO:LABELS DEFINITION
labels = zeros(height, width,3);
segmentation_labels = labels(:,:,1);
%O for bk 1 for fg
displacement_labels = labels(:,:,2:3);

%DATA TERM
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
%Bkg region is 3 vectors : RGB!
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

function result = SSDG (curr_patch, displaced_patch, gaussian)


end