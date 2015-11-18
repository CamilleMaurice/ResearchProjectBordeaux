function launchOneImage()
my_im=imread('Dataset/singleImage/planeSmall.jpeg');
roi=[2 2 2 2];
getSegmentation(my_im,my_im,roi)
end