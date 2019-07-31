%Pre-procesisng code for BIPGELS, June 4th samples.
%Written by Edna MB.For 5050 
%Last modification: July 10th-2019

close all
clear
clc

folder=['C:\Users\Ben\Documents\BIPJELS-Edna\5050\'];
%-->
test=['Test2and3\Series026'];
path_1=strcat(folder,test);
cd (path_1) %Folder where tif images are located
Files=dir('*ch01.tif'); %Read all of the images that end with ch01.tif (Usually Clay channel) 
%-->Resolution at x, y and z
resxy=0.08;
resz=0.299;

Z=length(Files);
%--> Initial and final layers
z1=3;
zend=90;
Perimeter=0;
aux=1;


for i=z1:zend

img=imread(Files(i).name);
img_i(:,:,aux)=img(225:end,225:end);


%% Lutidine
%-->
n=5;
FiltIm = imfilter(img,ones(n)/n^2,'symmetric'); %Apply a mean box filter
% for k = 1:3
%     FiltIm = imfilter(FiltIm,ones(n)/n^2,'symmetric');
% end

img_contrast=imadjust(FiltIm);

for k = 1:5
    img_contrast = imfilter(img_contrast,ones(n)/n^2,'symmetric'); %contrast enhancement
end

img_contrast_i(:,:,aux)=img_contrast;

%Thresholding opening-closing by reconstruction

    se = strel('disk',15);
    Ie = imerode(img_contrast,se); %Erosion of image 
    Iobr = imreconstruct(Ie,img_contrast); %Used to identify high-intensity objects in the mask

    % Opening-closing by reconstruction
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    %Iobrcbr =imcomplement(Iobrcbr);
    
    %--> Thresholding value
    n2=1.0;
    bw=imbinarize((Iobrcbr),graythresh(Iobrcbr)*n2); 
    %bw =~imbinarize(Iobrcbr,'adaptive','ForegroundPolarity','bright','Sensitivity',0.6); % As you reduce the '0.8' value, only brighter values are considered.

    %-->
    BWfinal= imopen(bw, strel('diamond',1)); %Remove small pores
    %BWfinal = imclose(BWfinal, strel('disk',5));
    BWfinal= bwareaopen(BWfinal,20,8) ; %Remove small pores
    BWfinal= ~bwareaopen(~BWfinal,20,8); %Remove small pores  

    BWfinal_i(:,:,aux)=BWfinal;  

    %%perimeter identification
    BWoutline=bwperim(BWfinal_i(:,:,aux));% last number indicates layer that you want to see
    Perimeter=Perimeter+sum(sum(BWoutline)); %total perimeter
    %displaying the video of the binary images
    Segout=img_contrast;
    Segout(BWoutline) = 255; 
    Segout_i(:,:,aux)=Segout;

    %-->
    imwrite(BWfinal_i(1:end,1:end,aux),strcat('C:\Users\Ben\Documents\BIPJELS-Edna\5050\Test2and3\Series026\Binary\Bipgel',num2str(i),'.tiff'));
    aux=aux+1;

end

implay(Segout_i);

%-->calculation of parameters
AreaOpenedIm_i=BWfinal_i;
[x,y,z]=size(AreaOpenedIm_i);
volume(1)=x*y*z*resxy^2*resz; %total volume in cubic micrometers
volume_sur=sum(sum(sum(AreaOpenedIm_i)))*resxy^2*resz;
surfacearea=Perimeter*resxy*resz;
volume(2)=surfacearea/volume_sur;%specific surface area
volume(3)=volume_sur/volume(1);%volume fraction

cd 'C:\Users\Ben\Documents\BIPJELS-Edna\'
%-->write in a new excel file
xlswrite('Results_concentration.xlsx',volume,'5050','F11:H11');



