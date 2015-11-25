classdef IndexClass
   properties
      index
      maxDisp
   end
   methods
       function obj = IndexClass( maxD )
          obj.maxDisp = maxD;
          obj.index = createIndex(obj);
       end
        function index = createIndex( obj )  
            maxDisplacement = obj.maxDisp;
         %'creating index'
            maxDisplacement
            TotLabels = 2*(maxDisplacement*2+1) * (maxDisplacement*2+1);
            index = zeros(TotLabels,4);
            for i = 1:1:TotLabels
                index(i,1) = i;
                %0 :: background, 1 :: object
                if ( i < TotLabels/2 +1 )
                    index(i,2) = 0;
                else
                    index(i,2) = 1;
                end
            end
        dx = -maxDisplacement;
        for i = 1:2*maxDisplacement+1:TotLabels
        
            if ( i == TotLabels/2 +1  )
                dx = -maxDisplacement;
            end
            
            cpt = 0;
            for dy = -maxDisplacement:1:maxDisplacement
                index(i+cpt, 4) = dy;
                cpt =cpt +1;
            end
            
            index(i:i+2*maxDisplacement+1, 3) = dx;
            dx = dx +1;
        end
        index = index(1:end -1, :);
        end
%         
%        function label = getLabel (index, seg, dx, dy)
%         
%         i4 = find(index(:, 4) == dy);
%         i3 = find(index(:, 3) == dx);
%         i2 = find(index(:, 2) == seg);
%         label = intersect(i4, i3);
%         label = intersect(label, i2);    
%         end
   end
end



