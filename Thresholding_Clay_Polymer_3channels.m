%Code for 3 channel confocal images
%Created by Allen
%Last modification: July 26th, 2019

close all hidden
clear
clc
warning off

%-->
cd 'G:\Allen\Aleksandra\U of A Fluorescent methacrylate\Series016 reg 4' %Folder where images are located


Files1 = dir('*ch00.tif'); %Read all of the images that end with ch00.tif (Usually Bitumen channel)
Files2 = dir('*ch01.tif'); %Read all of the images that end with ch01.tif (Usually Clay channel) 
Files3 = dir('*ch02.tif'); %Read all of the images that end with ch02.tif (Usually Polymer channel)

Z=length(Files1);%Read total layers of the images

for i = 1:Z
    imch1 = imread(Files1(i).name); %read Bitumen image at z=i
    imch2 = imread(Files2(i).name); %read Clay image at z=i
    imch3 = imread(Files3(i).name); %read Polymer image at z=i
    
    %crops the image
    imch1=imch1(1:512,1:512);
    imch2=imch2(1:512,1:512);
    imch3=imch3(1:512,1:512);
    
%% Bitumen Treatment
    %-->
    n = 4;
    FiltIm = imfilter(imch1,ones(n)/n^2,'symmetric'); %Apply an average box filter, symmetric means that 'Input array values outside the bounds of the array are computed by mirror-reflecting the array across the array border.'
    for k = 1:2 %the filter is applied 2 times or more to blur the image and obtain smoother edges   
        FiltIm = imfilter(FiltIm,ones(n)/n^2,'symmetric');
    end
    %thresholding factor rf=0.18;
    %BW = imbinarize(FiltIm,graythresh(FiltIm)*rf); %Binarize Image.
    BW= im2bw(FiltIm,0.18); %Binearize image. As you increase the numerical value, only brighter values are considered.
    bw_bitumen(:,:,i)=bwareaopen(BW,2); %Remove small objects, less than 2 pixels
    
%% Clay Treatment

    %cleaning and filtering
    %--> Level of detail of filtering
    m = 1;
    imch22 = imfilter(imch2,ones(m)/m^2); % Averaging filter. Other filter that can be used is the median filter 
    a = 0.18;b=0.8;c =1;d=0.9;
    
    %Thresholding  value can change throughout z
%     if i >=36 && i <=51
%          c = 0.52;d=0.8;
%          bw_clay1(:,:,i)=im2bw(Iobrcbr,a);
%     elseif i ==31
%         a=0.08; c=0.54;
%     bw_clay1(:,:,i)=imbinarize(Iobrcbr,graythresh(Iobrcbr)*c);
%     end
    bw_clay1(:,:,i) = im2bw(imch22,a);%% as you increase the 'a' value, only the brightest pixels are considered
    %Two different ways to binearize
    bw_clay2(:,:,i) = 1;%imbinarize(imch22,'adaptive','ForegroundPolarity','bright','Sensitivity',b);%As you reduce the 'b' value, only brighter values are considered.
    %Combines both binearuzations
    bw_clay(:,:,i)= bw_clay1(:,:,i).*bw_clay2(:,:,i);
    
%% Polymer Treatment
   
    %cleaning and filtering
    imch33 = imfilter(imch3,ones(3)/3^2); % Averaging filter. Other filter that can be used is the median filter 
    imch33 = imfilter(imch33,ones(3)/3^2);
    
    se3 = strel('disk',5);
    imch33_erode = imerode(imch33,se3); %Erosion of image 
    imch33_reconstruct = imreconstruct(imch33_erode,imch33); %Used to identify high-intensity objects in the mask 
    % External edges. Thresholding after opening-closing by reconstruction
    imch33_dilate = imdilate(imch33_reconstruct,se3);
    imch333 = imreconstruct(imcomplement(imch33_dilate),imcomplement(imch33_reconstruct));
    imch333 = imcomplement(imch333);
    
    a = 0.08;b=0.9;c =0.8;d=0.9;
    if i >=36 && i <=51
         c = 0.52;d=0.8;
         bw_clay1(:,:,i)=im2bw(Iobrcbr,a);
    elseif i ==31
        a=0.08; c=0.54;
    bw_clay1(:,:,i)=imbinarize(Iobrcbr,graythresh(Iobrcbr)*c);
    end
    bw_polymer1(:,:,i) = im2bw(imch333,a);%% as you increase the 'a' value, only the brightest pixels are considered
    bw_polymer2(:,:,i) = imbinarize(imch333,'adaptive','ForegroundPolarity','bright','Sensitivity',b);%As you reduce the 'b' value, only brighter values are considered.
    bw_polymer3(:,:,i) = 1;%imbinarize(Iobrcbr3,adaptthresh(Iobrcbr3,d,'ForegroundPolarity','bright'));
    bw_polymer4(:,:,i) = 1;%imbinarize(Iobrcbr3,graythresh(Iobrcbr3)*c);
    bw_polymer(:,:,i)= bw_polymer1(:,:,i).*bw_polymer2(:,:,i);
end

%cd 'F:\Allen\Aleksandra\U of A Fluorescent methacrylate\'

[bitumensq]=Bitumen_squish(Z,bw_bitumen);

%remove small objects
bw_clay=bwareaopen(bw_clay,1000,26);
bw_clay=~bwareaopen(~bw_clay,1000,26);
 
cd 'G:\Allen\Aleksandra\U of A Fluorescent methacrylate\Series016 reg 4\Binary'

for i =1:Z
    
%     bw_clay(:,:,i) = im2double(imread(strcat('Clay_bw',num2str(i),'.tiff')));
%     clay_bitsq(:,:,i)= im2double(imread(strcat('Clay_Final',num2str(i),'.tiff')));
    bw_bitumen(:,:,i) = imread(strcat('Bitumen_bw',num2str(i),'.tiff'));
    bitumensq(:,:,i) = im2double(imread(strcat('Bitumen_Final',num2str(i),'.tiff')));
    bw_polymer(:,:,i) = im2double(imread(strcat('Polymer_bw',num2str(i),'.tiff')));

end

bitumensq = bitumensq(1:512,1:512,:);
bw_polymer = bw_polymer(1:512,1:512,:);

clay_final=bw_clay-bitumensq;
clay_final(clay_final<0)=0;
polymer_final = bw_polymer-clay_final-bitumensq;
polymer_final(polymer_final<0)=0;

polymer_final=bwareaopen(polymer_final,1000,26);
polymer_final=~bwareaopen(~polymer_final,1000,26);

for i=1:Z
% imwrite((bw_bitumen(:,:,i)),strcat('Bitumen_bw',num2str(i),'.tiff')); %write the binary bitumen image to file
imwrite((bitumensq(:,:,i)),strcat('Bitumen_Final',num2str(i),'.tiff')); %%write the binary bitumen squished image to file
imwrite((bw_clay(:,:,i)),strcat('Clay_bw',num2str(i),'.tiff')); %%write the binary clay image to file
imwrite((clay_final(:,:,i)),strcat('Clay_Final',num2str(i),'.tiff')); %%write the final binary clay image to file
imwrite((bw_polymer(:,:,i)),strcat('Polymer_bw',num2str(i),'.tiff')); %%write the binary polymer image to file
imwrite((polymer_final(:,:,i)),strcat('Polymer_Final',num2str(i),'.tiff')); %%write the final binary polymer image to file

%% Water
% water(:,:,i)=imcomplement(bw_polymer); %assume everything except bitumen,clay, polymer in the sample is water
% imwrite(water(:,:,i),strcat('Water',num2str(i),'.tiff'))
end

binaryPath = 'G:\Allen\Aleksandra\U of A Fluorescent methacrylate\Series016 reg 4\Binary\';
rawPath = 'G:\Allen\Aleksandra\U of A Fluorescent methacrylate\Series016 reg 4\';

for i=1:Z
%     raw_bitumen(:,:,i)=imread(strcat(rawPath,Files1(i).name));   
%     outline_bitumen(:,:,i)=bwperim(bw_bitumen(:,:,i));    
%     outlineImg_bitumen(:,:,i)=raw_bitumen(:,:,i);
%     outlineImg_bitumen(outline_bitumen)=255;
    
    raw_clay(:,:,i)=imread(strcat(rawPath,Files2(i).name));   
    %bw_clay(:,:,i)=imread(strcat(binaryPath,'Clay_bw',num2str(i),'.tiff'));
    outline_clay(:,:,i)=bwperim(bw_clay(:,:,i));
    outlineImg_clay(:,:,i)=raw_clay(1:512,1:512,i);
    outlineImg_clay(outline_clay)=255;
   
    raw_polymer(:,:,i)=imread(strcat(rawPath,Files3(i).name));   
    %bw_polymer(:,:,i)=imread(strcat(binaryPath,'Polymer_bw',num2str(i),'.tiff'));
    outline_polymer(:,:,i)=bwperim(bw_polymer(:,:,i));
    outlineImg_polymer(:,:,i)=raw_polymer(1:512,1:512,i);
    outlineImg_polymer(outline_polymer)=255;
end

% handle1 = implay(outlineImg_bitumen);1
% handle1.Parent.Position = [200 50 1050 900];

handle2 = implay(outlineImg_clay);
handle2.Parent.Position = [200 50 1100 1100];

handle3 = implay(outlineImg_polymer);
handle3.Parent.Position = [1300 50 1100 1100];
























