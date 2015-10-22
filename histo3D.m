%3D histogram
%bins = 8
function histo = histo3D(image,n_bins)
w_bins=256/n_bins
[w,h,z] = size(image);
histo=zeros(n_bins,n_bins,n_bins);
for i=1:w
    for j=1:h
        image(i,j,:);
%         image(i,j,1)/w_bins
%         floor(image(i,j,1)/w_bins)
%         bin_r=floor(image(i,j,1)/w_bins)+1;
%         bin_r=(image(i,j,1)-mod(image(i,j,1),w_bins))/n_bins;
        bin_r=floor(double(image(i,j,1))/double(w_bins))+1;
        bin_g=floor(double(image(i,j,2))/double(w_bins))+1;
        bin_b=floor(double(image(i,j,3))/double(w_bins))+1;
%         bin_g=floor(image(i,j,2)/w_bins);
%         bin_b=floor(image(i,j,3)/w_bins);
        histo(bin_r,bin_g,bin_b) = histo(bin_r,bin_g,bin_b) +1;
    end
end
end