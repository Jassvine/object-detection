%% Initialization and uploading image
clear,clc,close all

th = 0.01;     % threshold for image segmentation
n  = 19;       % number of seeds per dimension
sideMax = 600; % maximum side length in pixels

% Load input image
[filename,pathname] = uigetfile('*.*');
Im = importdata([pathname,filename]);

% Resize image
scale = sideMax/(max(size(Im(:,:,1))));
Im = imresize(Im,scale*size(Im(:,:,1)));

%% Detecting the sign with bounding box
Mask = DetectSign(Im,th,n,sideMax); %detect sign function
bbox = regionprops(Mask,'BoundingBox');
bbox = bbox.BoundingBox;
SignImg = imcrop(Im,bbox);

%% subplot output 
figure('Position',[153,137,1610,780])
subplot(1,2,1),imshow(Im),title('Input image')
subplot(1,2,2),imshow(Im)
hold on,rectangle('Position',bbox,'EdgeColor','g','LineWidth',3)
hold off,title('Processed image')

figure(),imshow(SignImg),title('Detected sign')