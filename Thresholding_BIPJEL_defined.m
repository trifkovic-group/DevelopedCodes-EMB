%Pre-procesisng code for BIPGELS, June 4th samples.
%Written by Edna MB
%Last modification: July 9th-2019
close all hidden
clear
clc

%Code used for bipgels images binarization. Well defined structures.
%Code written by: Edna MB
%Last modification: July 31st,2019

%-->
folder=['C:\Users\Ben\Documents\BIPJELS-Edna\Surface3d-July25th\']; %Root directory
%-->
test=[]; %modify specific .tiff 
path_1=strcat(folder,test);
cd (path_1) %Folder where tif images are located
Files=dir('*ch01.tif'); %Read all of the images that end with ch01.tif (Usually Clay channel) 
%--> resolution of the image
resxy=0.059;
resz=0.243;

Z=length(Files);
%-->
z1=7; %initial layer
zend=41;%final layer
Perimeter=0;
aux=1;

%Binarization of each layer
for i=z1:zend
img_i(:,:,i)=imread(Files(i).name); %read specific layer

%--> Filter for smoothing boundaries and removing noise
n=7;

FiltIm = imfilter(img_i(:,:,i),ones(n)/n^2,'symmetric'); %Apply a box filter
for k = 1:2
    FiltIm = imfilter(FiltIm,ones(n)/n^2,'symmetric');
end
FiltIm_i(:,:,i)=FiltIm;

%Binarization
%--> Thresholding factor
n2=1.2;
BW_i(:,:,i) = imbinarize(FiltIm,graythresh(FiltIm)*n2); %Binarize Image
AreaOpenedIm_i(:,:,aux) = bwareaopen(BW_i(:,:,i),20); %Remove small objects
AreaOpenedIm_i(:,:,aux) = ~bwareaopen(~AreaOpenedIm_i(:,:,aux),20); %Remove small pores

%Used for displaying purposes
%perimeter identification
BWoutline=bwperim(AreaOpenedIm_i(:,:,aux));% Perimeter of the binary image
Perimeter=Perimeter+sum(sum(BWoutline));
Segout=img_i(:,:,i);
Segout(BWoutline) = 255; 
Segout_i(:,:,aux)=Segout;

%write binary images to use in Avizo later on
imwrite(AreaOpenedIm_i(:,:,aux),strcat(path_1,'\Binary\Bipjel',num2str(i),'.tiff'));
imwrite(~AreaOpenedIm_i(:,:,aux),strcat(path_1,'\Binary\Phase_2_',num2str(i),'.tiff'));
    
aux=aux+1;
end

implay(Segout_i);

%Calculation of total volume,specific surface area and volume fraction 
[x,y,z]=size(AreaOpenedIm_i);
volume(1)=x*y*z*resxy^2*resz; %total volume in cubic micrometers
volume_sur=sum(sum(sum(AreaOpenedIm_i)))*resxy^2*resz;
surfacearea=Perimeter*resxy*resz;
volume(2)=surfacearea/volume_sur;%specific surface area
volume(3)=volume_sur/volume(1);%volume fraction
xlswrite('Results_concentration.xlsx',volume,'7228','F11:H11'); %writes results in an specific excel file, sheet and cells range



%% Use for 3D printing, remove objects that are not connected to the top or bottom layers
%-->
% cd 'C:\Users\Ben\Documents\BIPJELS-Edna' %folder where remove_mid function is located
% Print=remove_mid(AreaOpenedIm_i);
% 
% aux=1;
% for i=z1:zend
% imwrite(Print(:,:,aux),strcat('C:\Users\Ben\Documents\BIPJELS-Edna\Surface3d-July25th\Print\Bipjel',num2str(i),'.tiff'));
% aux=aux+1;
% end

% cd 'C:\Users\Ben\Documents\BIPJELS-Edna\'
% %-->


% %% Outlined Original image/Clay
% % L=5;
%  %BWoutline=bwperim(BW_fin(:,:,L));% last number indicates layer that you want to see
%  BWoutline=bwperim(AreaOpenedIm_i(:,:,L));% last number indicates layer that you want to see
%  Segout=img_i(:,:,L);
%  Segout(BWoutline) = 255; 
%  figure, imshow((Segout)), title('Outlined original image');