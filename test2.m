function labels = test2(filename)

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
%imshow(rgb2gray(im))

%normalized greyscale image
%im = imread('car1.jpg');
im = imread(filename);
%GET INPUT IMAGE FROM GUI
m = double(rgb2gray(im))/255;
% figure(1)
% imshow(m)

%Convert image matrix into vector
reshape([m], [], 1);

[height,width] = size(m);

N = height*width;

% construct graph
E = edges4connected(height,width);

%prior avec distribution de couleurs avec histogrammes
%data term distance a lhistogramme
% prendre probabilite dun pxiel de couleurs dans lhistogramme.
% groupe de 16/8 pour la couleur (8*8*8)

%assuming the object is black and the background is white
costObj = abs(m);
costBkg = abs(1-m); 
imagesc(costObj)

%data term
T = [reshape([costObj], [], 1) , reshape([costBkg], [], 1) ];
T = sparse(T);

%smoothness term
alpha = 0.8 ;
sigma = 0.5 ;

Intensity_pq = [m(E(:,1)), m(E(:,2))];
Intensity_pq_diff = abs(Intensity_pq(:,1)-Intensity_pq(:,2)).*abs(Intensity_pq(:,1)-Intensity_pq(:,2));
%Repasser les E(:,) e
[idx,idy] = ind2sub([height, width], E(:,1));
[idx2,idy2] = ind2sub([height, width], E(:,2));

distance = sqrt((idx-idx2).^2 + (idy-idy2).^2);

B = (alpha*(exp(-Intensity_pq_diff))/(2*sigma*sigma))./distance;

A = sparse(E(:,1),E(:,2),B);

disp('calculating maximum flow');

[flow,labels] = maxflow(A,T);
labels = reshape(labels,[height width]);

%  figure(2)
%   imshow(m);
%  hold on;
%   %imagesc(labels); title(['labels // Alpha =',num2str(alpha),'sigma =', num2str(sigma)]);
%  colormap(gray)
%  imcontour(uint8(labels))
