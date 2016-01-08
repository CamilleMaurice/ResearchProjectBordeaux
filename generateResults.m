function [] = generateResults (folder, fileNumber, mask)
    gtPath = strcat(folder, 'ground-truth/');
    resultPath = strcat(folder, 'result/');
    
    if ~exist(resultPath, 'dir')
        mkdir(resultPath)
    end
% Write image to the right folder and write the csv file
    imwrite(mask, strcat(resultPath, sprintf('mask_%d.jpg', fileNumber-4)));
    
% The considered metrics is the number of different pixels between the
% ground truth and the mask
    
    gtlist = dir(gtPath);
    gtName = gtlist(fileNumber).name;
    gt = imread(strcat([gtPath, gtName]));       
    gt = gt(:,:,1);
    [h, w, ~] = size(gt);
    fileID = fopen(strcat(resultPath, 'result.csv'),'a');
    
    % Filenumber, number of different pixels, percentage of different pixels
    %diff = abs(sum(double(mask(:)) - double(gt(:))));
    gt = gt/255;
    mask = mask/255;
    diff = sum(sum(abs(double(gt)-double(mask)))) ;
    fprintf(fileID, sprintf(' %d, %d \n', diff, (diff*100)/(h*w)));
    
% [filtered ~] = getLargestCc(logical(imageObj.outMask),[],1);
% imwrite(filtered, sprintf('mask_%d.jpg',imageObj.fileNumber));
    
end