function labels = test2(filename,roi)

%Optical flow test
%addpath(genpath('OpticalFlow'));
%prev = imread('car1.jpg');
%curr = imread('car2.jpg');
%tic
%[u, v] = optic_flow_brox(prev, curr);
%[vx vy] = getOpticalFlow_Brox(prev,curr);
%toc
%figure(1)
%imagesc(sqrt(u.^2 + v.^2))
%GET INPUT IMAGE FROM GUI
im = imread(filename);
%m = double(rgb2gray(im))/255;
m = im;
[height,width] = size(m);
N = height*width;

%GET ROI COORDINATES FROM GUI
roi_im = imcrop(im,roi);

%GET BACKGROUND CONCATENATED PIXELS
RoiULx = roi(1);
RoiULy = roi(2);
RoiWidth = roi(3);
RoiHeight = roi(4);
s = size(im);

BkgL = imcrop(im,[1 1 RoiULx s(1)]);
BkgU = imcrop(im,[RoiULx, 1, RoiWidth, RoiULy]);
BkgD = imcrop(im,[RoiULx, RoiULy + RoiHeight, RoiWidth, s(2)-(RoiULy+RoiHeight)]);
BkgR = imcrop(im,[RoiULx+RoiWidth, 1, s(1)-(RoiULx), s(2)]);
%Check why BkgR is working it does not seem logical to me!

%vector froé rgb matrices B = reshape(permute(A,[3 1 2]),[],3);
BkgL = reshape([BkgL], [], 3);
BkgU = reshape([BkgU], [], 3);
BkgD = reshape([BkgD], [], 3);
BkgR = reshape([BkgR], [], 3);

Bkg = [BkgL; BkgU; BkgD ; BkgR] ;

%size(Bkg)
%compute 3d histograms

histo = histo3D(im,16);

%compute probabilities
reshape([m], [], 1);
% construct graph
E = edges4connected(height,width);

%prior avec distribution de couleurs avec histogrammes
%data term distance a lhistogramme
% prendre probabilite dun pxiel de couleurs dans lhistogramme.
% groupe de 16/8 pour la couleur (8*8*8)                                                                                                                                                                                                                                                                                                                                                                                                            
%assuming the object is black and the background is white
%costObj

costObj = -log(probsObj);
costBkg = -log(1-probsBkg); 


%data term
T = [reshape([costObj], [], 1) , reshape([costBkg], [], 1) ];
T = sparse(T);

%smoothness term
alpha = 0.8 ;
sigma = 0.5 ;

Intensity_pq = [m(E(:,1)), m(E(:,2))];
Intensity_pq_diff = abs(Intensity_pq(:,1)-Intensity_pq(:,2)).*abs(Intensity_pq(:,1)-Intensity_pq(:,2));
[idx,idy] = ind2sub([height, width], E(:,1));
[idx2,idy2] = ind2sub([height, width], E(:,2));

distance = sqrt((idx-idx2).^2 + (idy-idy2).^2);

B = (alpha*(exp(-Intensity_pq_diff))/(2*sigma*sigma))./distance;

A = sparse(E(:,1),E(:,2),B);

disp('calculating maximum flow');

[flow,labels] = maxflow(A,T);
labels = reshape(labels,[height width]);


