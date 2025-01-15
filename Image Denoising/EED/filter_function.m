%  Gaussian smoothening of the image filter of width sigma

function [smth] = filter_function(image, sigma);


smask = fspecial('gaussian', ceil(3*sigma), sigma);
smth = filter2(smask, image, 'same');