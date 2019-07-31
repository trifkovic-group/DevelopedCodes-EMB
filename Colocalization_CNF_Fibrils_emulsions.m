%Code written by Edna MB
%Thresholding CNF Fibril dispersions 
%Last modification date: July 30th-2019
%Parameters calculated as explaned in: https://www.physiology.org/doi/full/10.1152/ajpcell.00462.2010

clc
clear
close all hidden
warning off

Files=dir('*ch02.tif'); %Read the locations and names of the images that end with ch02.tif (Channel 3)
Files2=dir('*ch01.tif');%Read the locations and names of the images that end with ch01.tif (Channel 2)
Files3=dir('*ch00.tif');%Read the locations and names of the images that end with ch00.tif (Channel 1)

nlayers=1:1:length(Files);%Number of layers (images in File) 

for i=1:nlayers(end) %to evaluate each layer (each image)
   %Read images
    imch1=imread(strcat(Files(nlayers(i)).name));  %read Green channel image at z=i
    imch2=imread(strcat(Files2(nlayers(i)).name)); %read Blue channel image at z=i
    imch3=imread(strcat(Files3(nlayers(i)).name)); %read Red channel image at z=i
    
    %% Processing of channel 1 images
    %  Cleaning and Filtering
    imagech1=adapthisteq(imch1); % Enhances the contrast of the grayscale image. Other techniques (https://www.mathworks.com/help/images/contrast-enhancement-techniques.html)
    %Mean filter kernel size
    n=6;
    imagech1=imfilter(imagech1,ones(n)/n^2,'symmetric'); % Averaging filter. Other filter that can be used is the median filter 
    %thresholding factor
    %-->
    tf=0.5;
    imagech1_bw(:,:,i) = imbinarize(imagech1,graythresh(imagech1)*tf); %thresholding
    %-->
    imagech1_bw(:,:,i)= bwareaopen(imagech1_bw(:,:,i),50);%removing of small objects, less than 50 pixels
    imagech1_bit=imagech1;%enhanced image
    
    %For displaying purposes
    BWoutline=bwperim(imagech1_bw(:,:,i));% Perimeter of the binary image
    Segout=imagech1;%enhanced image
    Segout(BWoutline)= 255; %outlined perimeter
    Segout_i(:,:,i)=Segout; %store for video


    %-----------------------------------------%-------------------------------------------------%
    %% Processing of channel 2 images
     % Cleaning and Filtering
    imagech2=adapthisteq(imch2); % Enhances the contrast of the grayscale image. Other techniques (https://www.mathworks.com/help/images/contrast-enhancement-techniques.html)
    %second mean filter 
    %-->
    n2=6;
    imagech2=imfilter(imagech2,ones(n2)/n2^2,'symmetric'); % Averaging filter. Other filter that can be used is the median filter 
  
    %--> Thresholding
    tf2=0.7;
    imagech2_bw(:,:,i) = imbinarize(imagech2,graythresh(imagech2)*tf2);
    %remove small objects
    %-->
    imagech2_bw(:,:,i)= bwareaopen(imagech2_bw(:,:,i),50);
    imagech2_i(:,:,i)=(imagech2);
    imagech2_bit=imagech2;
    
    %For displaying purposes
    BWoutline=bwperim(imagech2_bw(:,:,i));% last number indicates layer that you want to see
    Segout2=imagech2;%enhanced image
    Segout2(BWoutline)= 255; %outlined perimeter
    Segout2_i(:,:,i)=Segout2; %store for video
    
    %% Processing of channel 3 images- I assume bitumen
    %Cleaning and Filtering
    imagech3=adapthisteq(imch3); % Enhances the contrast of the grayscale image. Other techniques (https://www.mathworks.com/help/images/contrast-enhancement-techniques.html)
    %Mean filter for smoothing edges
    %-->
    n3=10;
    FiltIm = imfilter(imagech3,ones(n3)/n3^2,'symmetric'); %Apply a box filter
    for k = 1:2
    FiltIm = imfilter(FiltIm,ones(n3)/n3^2,'symmetric');
    end
 %Thresholding
    tf3=0.15;
    BW = im2bw(FiltIm,tf3); %Binarize Image
    %-->
    se1 = strel('disk',3);%Structuring element 
    BW = imclose(~BW,se1);%Perform closing to remove small black spots in the particle network
    %-->
    BW = bwareaopen(~BW,500); %Remove small pores
    %watershed to separate touching objects
    imagech3_bw(:,:,i)= BW;

    BWoutline=bwperim(BW);% Perimeter of binary image
    % Used for colocalization calculations
    Bw_T=zeros(size(BW));
    Bw_T(BWoutline)=1;
    Bw_T1=imdilate(Bw_T,strel('disk',1));
    BWoutline_i(:,:,i)=Bw_T;

    %For displaying purposes
    Segout3=imagech3;
    Segout3(BWoutline) = 255; 
    Segout3_i(:,:,i)=Segout3;

    %location of perimeter of the bitumen droplets
    index=find(BW==1);
     
    imagech1_bit(index)=0;
    imagech2_bit(index)=0;
    %Formatting
    imagech1_i_bit(:,:,i)=double(imagech1_bit);
    imagech2_i_bit(:,:,i)=double(imagech2_bit);
    
    imwrite(imagech1_bw(:,:,i)-imagech3_bw(:,:,i),strcat('Binary_ch1','z_',num2str(i),'.tiff')); %write the binary images as tiff
    imwrite(imagech2_bw(:,:,i)-imagech3_bw(:,:,i),strcat('Binary_ch2','z_',num2str(i),'.tiff')); %write the binary images as tiff
    imwrite(imagech3_bw(:,:,i),strcat('Binary_ch3','z_',num2str(i),'.tiff')); %write the binary images as tiff
end

%Pre-procesing for colocalziation calculations
imagech1_i_bit_bw=imagech1_bw-imagech3_bw; %removes bitumen from channel 1
imagech1_i_bit_bw(imagech1_i_bit_bw==-1)=0;

imagech2_i_bit_bw=imagech2_bw-imagech3_bw-imagech1_i_bit_bw;%removes channel 3 and 1 from channel 2
imagech2_i_bit_bw(imagech2_i_bit_bw<0)=0;

%% colocalization parameters

Img_comb=imagech1_bw(:,:,:)+imagech2_bw(:,:,:);
nume=0;
denom1=0;
denom2=0;
nume_1=0;
denom11=0;
denom22=0;


meanch1=double(mean2(imagech1_i_bit(index)));
meanch2=double(mean2(imagech2_i_bit(index)));

%Calculation of Pearson's correlation coefficient and Manders overlap
%coefficient
   for i=1:length(index)
    %%PCC
        nume=nume+((imagech1_i_bit(index(i))-meanch1)*(imagech2_i_bit(index(i))-meanch2));
        denom1=denom1+(imagech1_i_bit(index(i))-meanch1)^2;
        denom2=denom2+(imagech2_i_bit(index(i))-meanch2)^2;   
    %%MOC    
        nume_1=nume_1+(imagech1_i_bit(index(i))*imagech2_i_bit(index(i)));
        denom11=denom11+(imagech1_i_bit(index(i)))^2;
        denom22=denom22+(imagech2_i_bit(index(i)))^2;
   end
   
   display(' Pearsons corralation factor PCC and Manders overlap coeffcient MOC')
   PCC= nume/sqrt(denom1*denom2)
   MOC=nume_1/sqrt(denom11*denom22)
   
   
   Ri_Total=sum(sum(sum(imagech1_i_bit_bw)));
   Gi_Total=sum(sum(sum(imagech2_i_bit_bw)));
   
   ri_cap=Img_comb.*(imagech1_i_bit_bw);
   ri_col=length(find(ri_cap==2));
   %Fractional overlap red channel
   display('Fractional overlap for channel 1')
   M1=ri_col/Ri_Total
   
   %Fractional overlap green channel
   gi_cap=Img_comb.*(imagech2_i_bit_bw);
   gi_col=length(find(gi_cap==2));
   display('Fractional overlap for channel 2')
   M2=gi_col/Gi_Total

%% Intersection on surface

%Calculates the fraction of the bitumen surface that is shared with channel
%1 and channel 2 objects
Total=length(find(BWoutline_i==1));
Int_ch1=length(find(imagech1_i_bit_bw+Bw_T1==2));
display('Fraction of bitumen surface colocalized with channel 1 objects')
F_ch1_bit=Int_ch1/Total
display('Fraction of bitumen surface colocalized with channel 2 objects')
Int_ch2=length(find(imagech2_i_bit_bw+Bw_T1==2));
F_ch2_bit=Int_ch2/Total

%% Videos of binary images
implay(Segout_i)
implay(Segout2_i)
implay(Segout3_i)