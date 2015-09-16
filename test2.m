function flow = test2()

% TEST2 Demonstrates gradient-based image
%   min-cut partitioning.
%
%   (c) 2008 Michael Rubinstein, WDI R&D and IDC
%   $Revision: 140 $
%   $Date: 2008-09-15 15:35:01 -0700 (Mon, 15 Sep 2008) $
%

im = imread('cow.jpeg');
%im = imread('Lenna.png');
figure(1)
imshow(rgb2gray(im))

%normalized greyscale image

m = double(rgb2gray(im)/255);

%Convert image matrix into vector
reshape([m], [], 1);

[height,width] = size(m);

disp('building graph');
N = height*width;

% construct graph
E = edges4connected(height,width);

%assuming the object is black and the background is white
costObj = abs(m);
costBkg = abs(1-m); 
imagesc(costObj)
%data term
T = [reshape([costObj], [], 1) , reshape([costBkg], [], 1) ];
T = sparse(T);

%smoothness term
alpha = 0.5 ;
sigma = 0.5 ;

Intensity_pq = [m(E(:,1)), m(E(:,2))];
Intensity_pq_diff = abs(Intensity_pq(:,1)-Intensity_pq(:,2)).*abs(Intensity_pq(:,1)-Intensity_pq(:,2));
distance = abs(E(:,1)-E(:,2));

%needs to be checked it doesnt change anything
B = (alpha*(exp(-Intensity_pq_diff))/(2*sigma*sigma))./distance;

A = sparse(E(:,1),E(:,2),B);

disp('calculating maximum flow');

[flow,labels] = maxflow(A,T);
labels = reshape(labels,[height width]);

 figure(2)
 imagesc(labels); title('labels // Alpha = 0.2 sigma = 1');
  %colormap(gray)
