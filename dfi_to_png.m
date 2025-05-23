clc;
clear all;
close all;

% Step 1: Load the background image
image_dfi = dfi2mat('scale_avg.dfi');
image = image_dfi.image;

% Step 2: Define the file name and path for the .png file
fileName = 'scale_avg.png';

% Step 3: Display the background image (optional)
figure;
imshow(image);
title('filename');

% Step 4: Save the background image as a .png file
imwrite(image, fileName);