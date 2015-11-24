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
index = createIndex();
labelCost = createLabelCost(index);
class = createClass(currentFrame, roi, index);

%APPEARANCE MODEL - UNARY/DATA TERM
histoBkg = histo3D( bkgRegion, n_bins );
histoObj = histo3D( reshape( objectRegion, [], 3 ), n_bins);

%Compute Probabilities - normalize histograms
probsObj = histoObj/size(reshape( objectRegion, [], 3 ),1);
probsBkg = histoBkg/size(bkgRegion,1);
% 
disp ('begin appearance model')
tic
Unary = zeros ( nLabels,nPixels );
UnaryMatrix = zeros( height, width, nLabels );
for label = 1:nLabels
    tic
    label
     
%      score1 = getApperanceSimilarity( label, index, windowOmega, currentFrame, previousFrame, maxDisplacement );
%      score2 = getApperanceModel( label, index, currentFrame, probsObj, probsBkg, n_bins);
   
   UnaryMatrix(:,:,label) = 0;%score1 + score2 ;   
   Unary(label, :) = reshape (UnaryMatrix(:,:,label), [], 1)';
   toc
end
toc
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

function index = createIndex()
'creating index'
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

function [Window] =  getNeigborhoodWindow ( y, x,image, WindowSize, maxDisplacement )

    half = floor(WindowSize/2);
    padSz = half + maxDisplacement;
    PaddedImg = padarray(image,[padSz, padSz],0);
    x_pad = x + padSz ;
    y_pad = y + padSz ;
    %do padarray to co ntrol borders
    Window = PaddedImg(y_pad - half : y_pad + half, x_pad - half : x_pad + half);
end

%Distance is defined as the number of different informations between 2
%labels. Max distance = 3; Min distance = 0 (if lp == lq)
function distance = getDistanceBtwLabels ( lp, lq, index )
    lp_info = index(lp,2:end);
    lq_info = index(lq,2:end);
    distance = sum(abs(lp_info - lq_info)); 
end

function score = getApperanceSimilarity( label, index, windowOmega, currentFrame, previousFrame, maxDisplacement )
%score will actually be the matrix of scores   

%BEWARE VARIABLES VISIBILITY
[height, width, ~] = size (currentFrame);

label_info = index(label, 2:end);
score = zeros(height, width);
dx = label_info(2);
dy = label_info(3);
gaussian = fspecial ( 'gaussian', [windowOmega, windowOmega], 0.5);

for y = 1:1:height
    for x = 1:1:width
        
        winCurr = getNeigborhoodWindow(y+dy, x+dx, currentFrame, windowOmega, maxDisplacement);
        winPrev = getNeigborhoodWindow(y, x, previousFrame, windowOmega, maxDisplacement);
        score(y,x) = GSAD ( winPrev, winCurr, gaussian );
    end
    if mod(y,10) == 0
        y
    end
end
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

function result = GSAD (currPatch, displacedPatch, gaussian)
%TO DO check if size of the curr_patch and displaced_patch is the same
[h,w,~] = size(currPatch);
img = abs(currPatch - displacedPatch);
result = sum(sum(filter2(gaussian, img)));
% for j = 1:h
%     for i = 1:w
%         result = result + gaussian(j,i) * abs(currPatch(j,i) - displacedPatch(j,i));
%     end
% end        

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


