function labels = getSegmentation(filename,roi)
%Returns the segmented image from the current frame and the initial
%selection of the object by the user : roi.

addpath(genpath('maxflow-v3.0'));
addpath(genpath('OpticalFlow'));

current_frame = imread( filename );
object_region = imcrop( current_frame, roi );
bkg_region = getBkgRegion( roi, current_frame );

[height,width,channels] = size( current_frame );
number_of_pixels = height*width;

%%LABELS DEFINITION
labels = zeros(height, width,2);
segmentation_labels = labels(:,:,1);
%O for bk 1 for fg
displacement_labels = labels(:,:,2);

%DATA TERM
%APPEARANCE SIMILARITY AS A FUNCTION OF the displacement
max_displacement = 2;
AS = zeros ( height, width, (2*max_displacement+1)^2 );
window_omega = 6;
std = 0.5
gaussian = fspecial('gaussian', window_omega, std);
cpt = 1;
for y = -max_displacement:max_displacement
    for x = -max_displacement:max_displacement
        AS(:,:,cpt) = getApperanceSimilarity( patch, displaced_patch );
        cpt = cpt +1;
    end
    
end
%Need a function to create window around each pixel of the image.


%%APPEARANCE MODEL
%Compute Histograms
%Bkg region is 3 vectors : RGB!
n_bins = 16;
histoBkg = histo3D( bkg_region, n_bins );
histoObj = histo3D( reshape( [object_region], [], 3 ), n_bins);

%Compute Probabilities - normalize histograms
probsObj = histoObj/size(reshape( [object_region], [], 3 ),1);
probsBkg = histoBkg/size(bkg_region,1);
%sum(sum(sum(probsObj))) == 1
%sum(sum(sum(probsBkg))) == 1

U_fg = zeros(height, width, 1);
U_bg = zeros(height, width, 1);
%OPTIMIZE THIS
% Rcolor = current_frame(:,:,1);
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


%SMOOTHNESS TERM
%MOTION COHERENCE
%ATTRIBUTE COHERENCE

reshape( [current_frame], [], 1 );
% construct graph
E = edges4connected( height, width );

%prior avec distribution de couleurs avec histogrammes
%data term distance a lhistogramme
% prendre probabilite dun pxiel de couleurs dans lhistogramme.
% groupe de 16/8 pour la couleur (8*8*8)

                                                                                                                                                                                                                                                                                                                                                                                                            
%assuming the object is black and the background is white

costObj = abs( 1 - current_frame );%;-log(probsObj);
costBkg = abs( current_frame );%-log(1-probsBkg); 


%data term
T = [reshape([costObj], [], 1) , reshape([costBkg], [], 1) ];
T = sparse(T);

%smoothness term
alpha = 0.8 ;
sigma = 0.5 ;

Intensity_pq = [current_frame(E(:,1)), current_frame(E(:,2))];
Intensity_pq_diff = abs(Intensity_pq(:,1)-Intensity_pq(:,2)).*abs(Intensity_pq(:,1)-Intensity_pq(:,2));
[idx,idy] = ind2sub([height, width], E(:,1));
[idx2,idy2] = ind2sub([height, width], E(:,2));

distance = sqrt((idx-idx2).^2 + (idy-idy2).^2);

B = (alpha*(exp(-Intensity_pq_diff))/(2*sigma*sigma))./distance;

A = sparse(E(:,1),E(:,2),B);

disp('calculating maximum flow');

[flow,labels] = maxflow(A,T);
labels = reshape(labels,[height width]);

end
function concatenatedImage = getBkgRegion(roi, image)
%This function has been unary tested
%Size of the image = size of the backgorund + size of the object

%GET BACKGROUND CONCATENATED PIXELS
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