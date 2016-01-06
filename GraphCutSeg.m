function [labels, binarySegmentationMask] = getSegmentationC1(currentFrame, previousFrame, roi)
addpath(genpath('GCMex'));

inputImage = imread('cow.jpeg');
%inputImage = imread('loutre.png');
inputImage = imresize(inputImage, 0.5);
figure(1)
imshow(inputImage)
inputImageVector = reshape(inputImage,[],3);
%user scribbles input as mask
inputBkg = imread('cowBkg.jpeg');
%inputBkg = imread('loutre_bg.png');
inputObj = imread('cowObj.jpeg');
%inputObj = imread('loutre_obj.png');

[height, width, ~] = size(inputImage);
numberPixels = height*width;

%Build the input matrices for maxflow
labelSetSize = 2;

%Class is the initialization
%Initialize to all pixels are Bkg.
class = zeros(1,numberPixels);
  class = randi(10,1, numberPixels);
  class(class <5) = 0;
  class(class >=5) = 1;

%unary = zeros(labelSetSize, numberPixels);

disp('Creation of the adjacency matrix and computation of pairwise costs')
tic 
%we set weigth where pixels are connected
adjMat = getAdjacencyMatrix(inputImage);
%pairWise = sparse(numberPixels, numberPixels);

PairWise = getAdjacencyMatrix(inputImage);
%p and q are neighboring pixels
[p,q,~] = find(adjMat);
for i = 1:numel(p)
        [neighbor1y neighbor1x] = ind2sub(size(inputImage),p(i));
        [neighbor2y neighbor2x] = ind2sub(size(inputImage),q(i));
        neighbor1 = inputImage(neighbor1y, neighbor1x);
        neighbor2 = inputImage(neighbor2y, neighbor2x); 
        PairWise(p(i),q(i)) = setPairWiseCost(neighbor1, neighbor2, neighbor1x,neighbor1y,neighbor2x,neighbor2y);
end
toc

%Set the labelCost
labelCost = [0 , 1 ; 1 , 0];

%Define unary potentials
%Find probabilities
n_bins = 16;
histoBkg = histo3D( reshape( inputBkg, [], 3 ), n_bins );
histoObj = histo3D( reshape( inputObj, [], 3 ), n_bins);

%Compute Probabilities - normalize
probsObj = histoObj/size(reshape( inputObj, [], 3 ), 1);
probsBkg = histoBkg/size(reshape( inputBkg, [], 3 ), 1);


unaryBkg = arrayfun(@(x1,x2,x3) setUnaryCost(x1,x2,x3,probsBkg), inputImageVector(:,1),inputImageVector(:,2), inputImageVector(:,3));
unaryObj = arrayfun(@(x1,x2,x3) setUnaryCost(x1,x2,x3,probsObj), inputImageVector(:,1),inputImageVector(:,2), inputImageVector(:,3));

%
figure(50)
imshow(255*reshape(unaryObj,[height, width]))

unaryBkg = unaryBkg';
unaryObj = unaryObj';

unary = [unaryBkg ; unaryObj];

%Check symmetry
%issym = @(x) isequal(x,x.')
%issym(PairWise)
% disp('Size of unary and PairWise')
% size(unary)
% size(PairWise)
disp('GCMex start')
[LABELS ENERGY ENERGYAFTER] = GCMex(class, single(unary), PairWise, single(labelCost), 1);
ENERGY
ENERGYAFTER
figure(2)
LABELS = reshape(LABELS, height, width);
imshow(LABELS)
%Create adjacency matrix
%http://stackoverflow.com/questions/3277541/construct-adjacency-matrix-in-matlab
end


%prob will be either probBkg or probObj
function cost = setUnaryCost (p1, p2, p3, prob)
n_bins = 16;
lambda = 20;
Rcolor = ceil(p1/(n_bins+1)) +1;
Gcolor = ceil(p2/(n_bins+1)) +1;
Bcolor = ceil(p3/(n_bins+1)) +1;
 
cost = -log( prob(Rcolor, Gcolor, Bcolor) +eps)*lambda;
%cost = prob(Rcolor, Gcolor, Bcolor)*1;
end


function cost = setPairWiseCost (p, q, xp, yp, xq, yq)
    intensityDiff = (p-q).^2;
    sigma = 1;
    %Equals to one in the case of the 4-neighborhood
    invDistance = 1/sqrt((xp - xq)^2 + (yp - yq)^2);
    if(invDistance ~= 1)
        disp('4-neighborhood, distance != 1 => BUG');
    end
    intensityDiff = double(intensityDiff);
    cost = exp(-intensityDiff/(2*sigma*sigma)) * 1;
end

function cost = setLabelCost(lp, lq)
    if ( lp == lq )
        cost = 0;
    else
        cost = 1;
    end

end

function W = getAdjacencyMatrix(I)

[m, n, ~] = size(I);

I_size = m*n;

% 1-off diagonal elements
V = repmat([ones(m-1,1); 0],n, 1);
V = V(1:end-1); % remove last zero

% n-off diagonal elements
U = ones(m*(n-1), 1);

% get the upper triangular part of the matrix
W = sparse(1:(I_size-1),    2:I_size, V, I_size, I_size)...
  + sparse(1:(I_size-m),(m+1):I_size, U, I_size, I_size);

% finally make W symmetric
W = W + W';
end