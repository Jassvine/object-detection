function [Part,Im] = DetectSign(Im,th,n,sideMax)
% This function takes in the  image Im
% and the other parameters then returns a 4-D
% matrix holding the parts

% % Parameters
% th = 0.01;    % threshold for image segmentation
% n  = 15;      % number of seeds per dimension
% sideMax = 600; % maximum side length in pixels 

% turn image to gray scale
if length(size(Im))>2
    Im = rgb2gray(Im);
end
% set the image size to suitable value
scale = sideMax/(max(size(Im)));
Im = imresize(Im,scale*size(Im));

% Create seeds grid
%basically subdividing input image and store in vector position
%and each pixel divided into parts
x = linspace(1,size(Im,2),n+2); %the length of the vector defined 
y = linspace(1,size(Im,1),n+2); %linspace func calculates interval to fit length
x = round(x); y = round(y); 

% seeding mask
MASK = false(size(Im));

% seeds
for i = 2:n+1
    for j = 2:n+1
        MASK(y(j),x(i)) = true;
    end
end

% Image proccesing
[y,x] = find(MASK);
Parts(:,:,1,1) = false(size(Im));   % preallocation
partsNum = 0; convexity = 0;        % preallocation

% start breaking down process
while ~isempty(y)
    mask = false(size(Im));
    mask(y(1),x(1)) = true;
    MASK(y(1),x(1)) = false;

    % Compute the weight array based on grayscale intensity differences.
    W = graydiffweight(Im, mask);

    % Segment the image using the weights.
    BW = imsegfmm(W, mask, th); %returns segmented  image
    % w is weights based on earlier W defined
    % mask is an array that defines seed locations 
    % th is treshold level defined earlier = 0.01

    % segment the part out of picture
    MASK(BW) = false;

    % Detect incorrect detected shapes
    BW = bwareafilt(BW,1); %bw extracts objects from binary image
    %the area of object is specified ina  range = 1
    CH = bwconvhull(BW); %generate convex hull of all objects in BW and generates CH 
     %CH  is a binary convex hull image
     %convex hull is basically the smallest convex polygon containing all
     %the given points
     %CH is important to manage  flexibility in capturing shape of an object in an image
    partsNum = partsNum+1;
    Peri = regionprops(BW,'Perimeter','Area'); %returns measurements 
    if Peri.Area > (numel(BW)/100) %numel is to return number of elements in array BW
        Peri = Peri.Perimeter;
        PeriConv = regionprops(CH,'Perimeter');
        PeriConv = PeriConv.Perimeter;
        convexity(partsNum) = PeriConv/Peri;
    else
        convexity(partsNum) = 0;
    end

    % save it in parts
    Parts(:,:,partsNum) = BW;
    [y,x] = find(MASK);
end

% Selecting part with highest convexity
[~,idx] = max(convexity);
Part = Parts(:,:,idx);



