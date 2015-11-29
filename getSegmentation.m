function [labels, binarySegmentationMask] = getSegmentation(currentFrame, previousFrame, roi)
lambda1=0.2;
lambda2=0.2;
lambda3=1;
Ec=2;
n_bins = 16;
maxDisplacement = 2;%in the paper they used 16
windowOmega = 7;
%should be renamed
std = 0.5;
nLabels = 2*(2*maxDisplacement+1)*(2*maxDisplacement+1);
gaussian = fspecial('gaussian', windowOmega, std);

% addpath(genpath('maxflow-v3.0'));

% currentFrame = imread( currFilename );
% previousFrame = imread( prevFilename );
imageObj = ImageClass(currentFrame, roi);

% objectRegion = imcrop( currentFrame, roi );
% bkgRegion = getBkgRegion( roi, currentFrame );

[height, width, ~] = size(currentFrame);
nPixels = height*width;

%LABELS DEFINITION
%index = createIndex();

indexObj = IndexClass(maxDisplacement);
index = indexObj.index;
nLabels = size(index, 1);
labelCost = createLabelCost(indexObj);
class = createClass(imageObj, indexObj);

%APPEARANCE MODEL - UNARY/DATA TERM

%Compute Probabilities - normalize histograms
probsObj = getProbabilities(imageObj,'fg',n_bins);
probsBkg = getProbabilities(imageObj,'bg',n_bins);

disp ('begin appearance model')
Unary = zeros ( nLabels,nPixels );
UnaryMatrix = zeros( height, width, nLabels );
for label = 1:nLabels
    
   score1 = getApperanceSimilarity( label, indexObj, imageObj, gaussian, previousFrame );
   score2 = getApperanceModel( label, index, currentFrame, probsObj, probsBkg, n_bins);

   UnaryMatrix(:,:,label) = double(score1) + score2 ;   
   Unary(label, :) = reshape (UnaryMatrix(:,:,label), [], 1)';
   
end
disp (' appearance model done')
%APPEARANCE MODEL - UNARY/DATA TERM DONE

%CHECK: shouldn't unary be of size nLabels*nPixels? 
%It must be nLabels*nbPixels :)

%DONE:SMOOTHNESS TERM spatial
Spatial_Pairwise=zeros(nLabels,nLabels);

%DONE:MOTION COHERENCE
%for now we use SAD to compute distance 
DxDy=index(1:nLabels/2,3:4);
D=abs(DxDy(:,1)-DxDy(:,2));
D_LpLq=zeros(nLabels/2,nLabels/2);
for i=1:length(D)
    K=D(i)+D;
    D_LpLq(i,:)=K';
end
D_LpLq=lambda3*D_LpLq;
Spatial_Pairwise = [D_LpLq D_LpLq; D_LpLq D_LpLq];

%DONE:ATTRIBUTE COHERENCE
Spatial_Pairwise(nLabels/2+1:end,1:nLabels/2) = Spatial_Pairwise(nLabels/2+1:end,1:nLabels/2)+Ec;
Spatial_Pairwise(1:nLabels/2,nLabels/2+1:end) = Spatial_Pairwise(1:nLabels/2,nLabels/2+1:end)+Ec;
%fitting the smooth term init
Spatial_Pairwise=reshape(Spatial_Pairwise,[],1);

%SPATIAL SMOOTHNESS TERM DONE

reshape( currentFrame, [], 1 );

%TODO TUNE GCMEX FOR SMOOTHNESS TERM
%TODO CALL GC MEX
%gfet the labels
%TODO: construct the binary segmentation mask from the output labels
% nPixels
% size(class)
% size(Unary)
% size(Spatial_Pairwise)
% size(labelCost)
[labels, energy, energyafter] = GCMex(class', single(Unary), Spatial_Pairwise, single(labelCost),0,width,height);
% labels=zeros(height,width);
% labels( 10:20,10:20)=1;
%   display('the energy is')
  energy
   energyafter
    binarySegmentationMask=fromLabelstoMask(labels, width, height);
    roi =  fromMaskToRoi ( labels, width );
end

% For initialization all dx, dy are set to 0, we choose the label
% corresponding to bg/fg with displacement 00;
function class = createClass (imageObj, indexObj)
    imagetest = imageObj.image;
    [h, w, ~] = size(imagetest);
    class = zeros (h, w);
    labelFg00 = getLabel(indexObj, 1, 0, 0);
    labelBg00 = getLabel(indexObj, 0, 0, 0);
    
    for y = 1:h
        for x = 1:w
            if (isInRoi(imageObj, y, x))
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

%Create a label cost matrix according to the format wanted by GCMEX.
function labelCost = createLabelCost ( indexObj )
    index = indexObj.index;
    nLabels = size (index, 1); 
    
    labelCost = zeros (nLabels, nLabels);
    for lp = 1:nLabels
        for lq = 1:nLabels
            labelCost(lp,lq) = getDistanceBtwLabels( indexObj, lp, lq );
        end
    end    
end

function score = getApperanceSimilarity( label, indexObj, imageObj, gaussian, previousFrame )
%score will actually be the matrix of scores   

index = indexObj.index;
maxDisplacement = indexObj.maxDisp;
currentFrame = imageObj.image;

[height, width, ~] = size (currentFrame);
score = zeros(height, width, 3);

label_info = index(label, 2:end);
dx = label_info(2);
dy = label_info(3);

currentFrameDisplaced = displaceImage( imageObj, maxDisplacement, dy, dx );

score(:,:,1) = abs(conv2 (double(currentFrameDisplaced(:,:,1)), double(gaussian), 'same') - conv2(double(previousFrame(:,:,1)),double (gaussian), 'same'));
score(:,:,2) = abs(conv2 (double(currentFrameDisplaced(:,:,2)), double(gaussian), 'same') - conv2(double(previousFrame(:,:,2)),double (gaussian), 'same'));
score(:,:,3) = abs(conv2 (double(currentFrameDisplaced(:,:,3)), double(gaussian), 'same') - conv2(double(previousFrame(:,:,3)),double (gaussian), 'same'));

score = score(:,:,1) + score (:,:,2) + score(:,:,3);
end

function U_score = getApperanceModel( label, index, currentFrame, probsObj, probsBkg,n_bins )
    [height, width, ~] = size(currentFrame);
    
    if ( sum(ismember(find(index(:,2)==0),label) ) )
        isBkg = true(1);
    else
        isBkg = false(1);
    end
    
    U_score = zeros(height, width, 1);
    % TO DO : OPTIMIZE THIS BY THE CALL OF ANONYM FUNCTION OVER A MATRIX
    for j = 1 : height
        for i = 1 : width
            %OPTIMIZE THIS
            Rcolor = floor(currentFrame(j,i,1)/n_bins)+1;
            Gcolor = floor(currentFrame(j,i,2)/n_bins)+1;
            Bcolor = floor(currentFrame(j,i,3)/n_bins)+1;
            
            %If the label is BG
            if ( isBkg )
                U_score(j,i) = -log (probsBkg(Rcolor,Gcolor,Bcolor) +eps);
            else
                U_score(j,i) = -log (probsObj(Rcolor,Gcolor,Bcolor) +eps); 
            end
        end
    end
end

function roi = fromMaskToRoi ( mask, width )
CC=bwconncomp(mask);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
%[col row]=ind2sub(width,height);
bb=regionprops(CC,'BoundingBox');%not good!! for several bb
roi=[bb.BoundingBox(1:2)+0.5 bb.BoundingBox(3:4)];

end

function mask = fromLabelstoMask(labels,height,width)
mask=zeros(height,width);
med=max(max(labels))/2;
mask(labels>med)=1;%TODO:CH:check the strict superiority

end
