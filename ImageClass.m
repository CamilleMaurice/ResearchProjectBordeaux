classdef ImageClass
    properties
        path
        fileNumber
        image
        rect
        fgVector
        bgVector
        mask
    end
    methods
          %Constructor
          function obj = ImageClass( folder, fileNumber, curr, roi )
          obj.path = folder;    
          obj.fileNumber = fileNumber;
          obj.image = curr;
          obj.rect = roi;
          obj.fgVector = getObjVector(obj);
          obj.bgVector = getBgVector(obj);
          end
          
          %Returns the object region in a vector
          function fg = getObjVector(obj)
             fg = imcrop( obj.image, obj.rect );   
             fg = reshape( fg, [], 3 );
          end
          
          %Returns the background region in a vector
          function bg = getBgVector(obj)
            %This function has been unary tested
            %Size of the image = size of the backgorund + size of the object
            size_im = size(obj.image);
            roi = obj.rect;
            RoiULx = roi(1); RoiULy = roi(2); RoiWidth = roi(3); RoiHeight = roi(4);
            
            BkgL = imcrop(obj.image,[1, 1, RoiULx-1, size_im(1)]);
            BkgU = imcrop(obj.image,[RoiULx, 1, RoiWidth-1, RoiULy-1]);
            BkgD = imcrop(obj.image,[RoiULx, RoiULy + RoiHeight + 1, RoiWidth-1, size_im(1)- (RoiULy + RoiHeight +1)]);
            BkgR = imcrop(obj.image,[RoiULx+RoiWidth+1, 1, size_im(2)-(RoiULx + RoiWidth +1), size_im(1)]);
            
            %from matrix to vector;
            BkgL = reshape(BkgL, [], 3);
            BkgU = reshape(BkgU, [], 3);
            BkgD = reshape(BkgD, [], 3);
            BkgR = reshape(BkgR, [], 3);
            
            bg = [BkgL; BkgU; BkgD ; BkgR] ;
          end
          
          %Returns probabilities
          function probs = getProbabilities (obj, type, nbins)
              if ( strcmp(type ,'bg') == 1 )
                    region = obj.bgVector;
              elseif (strcmp(type ,'fg') == 1 )
                      region = obj.fgVector;
              else
                  disp('Usage error, type must be either fg or bg')                  
              end
              histo = histo3D( region, nbins );
              probs = histo/size( region, 1 );              
          end
          
          %Returns true is the coordinate (y,x) is inside the ROI, else false.
          function  [bool]  = isInRoi (obj, y, x)
              
              roi = obj.rect;
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
          
          %Returns a displaced image along a given direction (dy, dx)
          function imgD = displaceImage ( obj, maxD, dy, dx )
              img = obj.image;
              paddedImg = padarray( img,[maxD, maxD], 'replicate' );
              imgD = paddedImg( maxD+1 +dy:end-maxD +dy, maxD+dx+1:end-maxD+dx, : );
          end
          
          %Update of the model, trigger generation of the results or view
          function [] = updateMask (obj, segmentation)
              obj.mask = segmentation;
              folder = obj.path;
              generateResults(folder, obj.fileNumber, obj.mask);              
          end
    end
end

