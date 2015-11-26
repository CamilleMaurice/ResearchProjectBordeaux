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
objectRegion = imcrop( currentFrame, roi );
bkgRegion = getBkgRegion( roi, currentFrame );

[height, width, ~] = size(currentFrame);
sizeIm = [height, width];
nPixels = height*width;

%LABELS DEFINITION
%index = createIndex();

indexObj = IndexClass(maxDisplacement);
index = indexObj.index;
%nLabels = size(index, 1)
labelCost = createLabelCost(indexObj);
class = createClass(currentFrame, roi, indexObj);

%APPEARANCE MODEL - UNARY/DATA TERM
histoBkg = histo3D( bkgRegion, n_bins );
histoObj = histo3D( reshape( objectRegion, [], 3 ), n_bins);

%Compute Probabilities - normalize histograms
probsObj = histoObj/size(reshape( objectRegion, [], 3 ),1);
probsBkg = histoBkg/size(bkgRegion,1);
% 
disp ('begin appearance model')

Unary = zeros ( nLabels,nPixels );
UnaryMatrix = zeros( height, width, nLabels );
for label = 1:nLabels
    
   score1 = getApperanceSimilarity( label, index, windowOmega, currentFrame, previousFrame, maxDisplacement );
   score2 = getApperanceModel( label, index, currentFrame, probsObj, probsBkg, n_bins);

   UnaryMatrix(:,:,label) = double(score1) + score2 ;   
   Unary(label, :) = reshape (UnaryMatrix(:,:,label), [], 1)';
   
end

disp (' appearance model done')
%APPEARANCE MODEL - UNARY/DATA TERM DONE

%CHECK: shouldn't unary be of size nLabels*nPixels?

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
BkgL = reshape(BkgL, [], 3);
BkgU = reshape(BkgU, [], 3);
BkgD = reshape(BkgD, [], 3);
BkgR = reshape(BkgR, [], 3);

concatenatedImage = [BkgL; BkgU; BkgD ; BkgR] ;
end

% For initialization all dx, dy are set to 0, we choose the label
% corresponding to bg/fg with displacement 00;
function class = createClass (image, roi, indexObj)
    [h, w, ~] = size(image);
    class = zeros (h, w);
    labelFg00 = getLabel(indexObj, 1, 0, 0);
    labelBg00 = getLabel(indexObj, 0, 0, 0);
    
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

%Create a displaced image along a given direction
function imgD = displaceImage ( img, maxD, dy, dx )
    paddedImg = padarray( img,[maxD, maxD], 'replicate' );
    imgD = paddedImg( maxD+1 +dy:end-maxD +dy, maxD+dx+1:end-maxD+dx, : ); 
end

%Compute the integralImage to fasten computation
function res = integralImage ( img )
    RChannel = img(:,:,1); GChannel = img(:,:,2); BChannel = img(:,:,3);
    resR = cumsum(cumsum(RChannel')');
    resG = cumsum(cumsum(GChannel')');
    resB = cumsum(cumsum(BChannel')');
    res = resR + resG + resB;
end

function [Window] =  getNeigborhoodWindow ( y, x, image, WindowSize, maxDisplacement )

    half = floor(WindowSize/2);
    padSz = half + maxDisplacement;
    PaddedImg = padarray(image,[padSz, padSz],0);
    x_pad = x + padSz ;
    y_pad = y + padSz ;
    %do padarray to co ntrol borders
    Window = PaddedImg(y_pad - half : y_pad + half, x_pad - half : x_pad + half);
end


function score = getApperanceSimilarity( label, index, windowOmega, currentFrame, previousFrame, maxDisplacement )
%score will actually be the matrix of scores   

%BEWARE VARIABLES VISIBILITY
[height, width, ~] = size (currentFrame);

label_info = index(label, 2:end);
score = zeros(height, width);
dx = label_info(2);
dy = label_info(3);

currentFrameDisplaced = displaceImage( currentFrame, maxDisplacement, dy, dx);
currentFrameDisplacedIntegral = integralImage( double (currentFrameDisplaced) );
previousFrameIntegral = integralImage( double(previousFrame) );
score = abs(currentFrameDisplaced - previousFrame);
score = score(:,:,1) + score (:,:,2) + score(:,:,3);
%gaussian = fspecial ( 'gaussian', [windowOmega, windowOmega], 0.5);
%  for y = 1:windowOmega:height
%      for x = 1:windowOmega:width   
%           
%           winCurr = getNeigborhoodWindow(y, x, abs(currentFrameDisplacedIntegral-previousFrameIntegral), windowOmega, maxDisplacement);
%           score(y:y+windowOmega,x:x+windowOmega) = winCurr(1,1) + winCurr(windowOmega, windowOmega) - winCurr(1, windowOmega) - winCurr(windowOmega, 1);
%             
%      end   
%  end


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
                U_score(j,i) = -log (probsBkg(Rcolor,Gcolor,Bcolor));
            else
                U_score(j,i) = -log (probsObj(Rcolor,Gcolor,Bcolor)); 
            end
        end
    end
end


% function result = GSAD (currPatch, displacedPatch, gaussian)
% %TO DO check if size of the curr_patch and displaced_patch is the same
% [h,w,~] = size(currPatch);
% img = abs(currPatch - displacedPatch);
% result = sum(sum(filter2(gaussian, img)));
% % for j = 1:h
% %     for i = 1:w
% %         result = result + gaussian(j,i) * abs(currPatch(j,i) - displacedPatch(j,i));
% %     end
% % end        
% 
% end

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


