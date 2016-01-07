function [ ] = getSegmentationC1( imageObj )
addpath(genpath('GCMex'));

inputImage = imageObj.image;
inputImage = imresize(inputImage, 0.5);

inputImageVector = reshape(inputImage,[],3);
[height, width, ~] = size(inputImage);
numberPixels = height*width;

%Class is the initialization
%Initialize to all pixels are Bkg.
class = zeros(1,numberPixels);

 disp('Creation of the adjacency matrix and computation of pairwise costs')
 tic 
%we set weigth where pixels are connected
PairWise = getAdjacencyMatrix(inputImage);
%p and q are neighboring pixels
[p,q,~] = find(PairWise);
for i = 1:numel(p)
        [neighbor1y neighbor1x] = ind2sub(size(inputImage),p(i));
        [neighbor2y neighbor2x] = ind2sub(size(inputImage),q(i));
        neighbor1 = inputImage(neighbor1y, neighbor1x);
        neighbor2 = inputImage(neighbor2y, neighbor2x); 
        PairWise(p(i),q(i)) = setPairWiseCost(neighbor1, neighbor2, neighbor1x,neighbor1y,neighbor2x,neighbor2y);
end
disp ('time to build pairwise')
toc

%Set the labelCost
labelCost = [0 , 1 ; 1 , 0];

%Define unary potentials

%Compute Probabilities
nbins = 16;
probsObj = imageObj.getProbabilities('fg', nbins);
probsBkg = imageObj.getProbabilities('bg', nbins);

unaryBkg = arrayfun(@(x1,x2,x3) setUnaryCost(x1,x2,x3,probsBkg), inputImageVector(:,1),inputImageVector(:,2), inputImageVector(:,3));
unaryObj = arrayfun(@(x1,x2,x3) setUnaryCost(x1,x2,x3,probsObj), inputImageVector(:,1),inputImageVector(:,2), inputImageVector(:,3));

unaryBkg = unaryBkg';
unaryObj = unaryObj';

unary = [unaryBkg ; unaryObj];

%Check symmetry
%issym = @(x) isequal(x,x.')
%issym(PairWise)

disp('GCMex start')
[labels ENERGY ENERGYAFTER] = GCMex(class, single(unary), PairWise, single(labelCost), 1);

% Update the model, the ImageClass object

segmentation = reshape(labels, height, width);
segmentation = imresize(segmentation, 2);
imageObj.updateMask(segmentation); 

    
end


%prob will be either probBkg or probObj
function cost = setUnaryCost (p1, p2, p3, prob)
n_bins = 16;
lambda = 1;
Rcolor = ceil(p1/(n_bins+1)) +1;
Gcolor = ceil(p2/(n_bins+1)) +1;
Bcolor = ceil(p3/(n_bins+1)) +1;
 
cost = -log( prob(Rcolor, Gcolor, Bcolor) +eps)*lambda;
%cost = prob(Rcolor, Gcolor, Bcolor)*1;
end


function cost = setPairWiseCost (p, q, xp, yp, xq, yq)
    intensityDiff = (p-q).^2;
    sigma = 200;
    %Equals to one in the case of the 4-neighborhood
    invDistance = 1/sqrt((xp - xq)^2 + (yp - yq)^2);
    if(invDistance ~= 1)
        disp('4-neighborhood, distance != 1 => BUG');
    end
    intensityDiff = double(intensityDiff);
    cost = exp(-intensityDiff/(2*sigma*sigma)) * 1;
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