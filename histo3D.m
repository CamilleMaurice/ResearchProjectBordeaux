%3D histogram
%bins = 8
function histo = histo3D(image,n_bins)
w_bins=256/n_bins;
[l,z] = size(image)
histo=zeros(n_bins,n_bins,n_bins);
for i=1:l
        bin_r=floor(double(image(l,1))/double(w_bins))+1;
        bin_g=floor(double(image(l,2))/double(w_bins))+1;
        bin_b=floor(double(image(l,3))/double(w_bins))+1;
        histo(bin_r,bin_g,bin_b) = histo(bin_r,bin_g,bin_b) +1;
end
end