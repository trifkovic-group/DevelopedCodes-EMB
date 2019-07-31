%Code written by Edna MB
%Thresholding of Bitumen/Clay/Polymer images taken with confocal microscopy
%Last modification date: February 14th-2019

close all
clear
clc

addpath('D:\Edna\Example with all codes') 
cd Series002 %Folder where tif images are located
Files=dir('*ch01.tif'); %Read all of the images that end with ch00.tif (Usually Clay channel) 
Files2=dir('*ch00.tif');%Read all of the images that end with ch01.tif (Usually Bitumen channel)

ak=1:1:length(Files);%Number of layers (images in File) 

for i=1:ak(end) %to evaluate each layer (each image)
   
    imch1=imread(strcat(Files(ak(i)).name));  %read Clay image at z=i
    imch2=imread(strcat(Files2(ak(i)).name)); %read Bitumen image at z=i
    
    %% Processing of Clay images
    % Cleaning and Filtering
    image22=adapthisteq(imch1); % Enhances the contrast of the grayscale image. Other techniques (https://www.mathworks.com/help/images/contrast-enhancement-techniques.html)
    image22=imfilter(image22,ones(3)/3^2); % Averaging filter. Other filter that can be used is the median filter 
    
    %-----------------------------------------%-------------------------------------------------%
    %Determine Voids
    
    se = strel('disk',8);
    Io = imopen(image22,se); %Opening of image with structure SE
    
    Ie = imerode(image22,se); %Erosion of image 
    Iobr = imreconstruct(Ie,image22); %Used to identify high-intensity objects in the mask
    Ioc = imclose(Io,se); %Opening closing
    % Opening-closing by reconstruction
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    %Threshholding
    bw = imbinarize(Iobrcbr,'adaptive','ForegroundPolarity','bright','Sensitivity',0.85); % As you reduce the '0.8' value, only brighter values are considered.
    
    %-----------------------------------------%-------------------------------------------------%
    %Determine Edges
    %Explanation cen be found in: https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
    
    [~, threshold] = edge(image22, 'sobel'); %Detects all the elements using a gradient function
    fudgeFactor = 0.3; %Trial and error
    BWs1 = edge(image22,'sobel', threshold* fudgeFactor);
    BWs=BWs1;
    
    % Dilate the image
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    BWsdil = imdilate(BWs,[se90 se0]);
    
    
    % Fill interior gaps
    BWdfill = bwareaopen(BWsdil, 100);
    BWdfill=imfill(BWdfill, 'holes');
    
    %Smoothing the object
    seD = strel('diamond',1);
    BWfinal = imerode(BWdfill,seD);
    
    %-----------------------------------------%-------------------------------------------------%
    
    bw=imreconstruct(BWfinal,bw);%Merge information of edges and voids
    bw_i(:,:,i)=bw; % Saves a images ina 3D array to be used in the Box-counting method 
    imwrite(bw,strcat('Clay',num2str(i),'.tiff')); %Writes binary images as .tiff
      
    %% Processing of bitumen images
    
%Normal Thresholding
    BWfinalgreen=im2bw(imch2,0.4); % as you increase the second value [0 1], only the brightest pixels are considered 
    BWfinalgreen=imfill(BWfinalgreen,'holes');
    BWfinalgreen=bwareaopen(BWfinalgreen,50,26);
    
   bw_bitumen(:,:,i)=BWfinalgreen;% SUPER IMPORTANT: This array is the one used in the squishing part 
   imwrite(bw_bitumen(:,:,i),strcat('Bitumen',num2str(i),'.tiff'));
    
end
    
cd 'D:\Edna\Example with all codes' %Go back to the previous folder where the function is located

[I3,Z]=Bitumen_squish_Fr(bw_bitumen);

 %% This section removes the squished Bitumen from the Clay images
 
for i=1:Z
   Clay=imread(strcat('Clay',num2str(i),'.tiff')); 
   Bitumen=imread(strcat('Slice_cut',num2str(i),'.tiff'));
   Clay_bit(:,:,i)=Clay-logical(Bitumen);
   Clay_bit_plus(:,:,i)=Clay+logical(Bitumen);
   imwrite(Clay_bit(:,:,i),strcat('Clay_bit',num2str(i),'.tiff')); 
   imwrite(Clay_bit_plus(:,:,i),strcat('Clay_bit_plus',num2str(i),'.tiff')); 
end

%% Calculate 3D-Fractal dimension

resxy=0.180; %Resolution xy (the one LASX gives you)
resz=0.122;  %Resolution z

Df_3d=Box_Counting_3D(resxy,resz,Z,Clay_bit);

%% Outlined Original image/Clay
% BWoutline=bwperim(bw_i(:,:,23));% last number indicates layer that you want to see
% Segout=imread(strcat(Files(ak(23)).name));
% Segout(BWoutline) = 255; 
% figure, imshow((Segout)), title('Outlined original image');