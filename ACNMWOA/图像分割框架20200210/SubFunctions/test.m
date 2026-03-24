clc;
clear;
close all;
img_14037=load('14037.mat');
img_500_14037=load('500_14037.mat');
a=img_14037.imageLabelCell{1, 1};
b=img_500_14037.groundTruth{1, 1}.Segmentation;
c=double(b);

figure;
X=label2rgb(a);
imshow(X);
figure;
Y=label2rgb(b);
imshow(Y);
figure;
Y=label2rgb(c);
imshow(Y);


