function histo = histo3D(image,n_bins)
w_bins = 256/n_bins;
[l,z] = size(image);
histo = zeros(n_bins,n_bins,n_bins);
for i=1:l
        bin_r = floor(double(image(i,1))/double(w_bins))+1;
        bin_g = floor(double(image(i,2))/double(w_bins))+1;
        bin_b = floor(double(image(i,3))/double(w_bins))+1;
        histo(bin_r,bin_g,bin_b) = histo(bin_r,bin_g,bin_b) +1;
end

%%Unary tested
%%size(image) == sum(sum(sum(histo)))
end