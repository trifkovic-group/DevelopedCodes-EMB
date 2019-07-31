close all
clear
clc

%-->
addpath('D:\Edna\Example with all codes') 
cd Series002

%% Read the 2D-Image 
%-->
bw=logical(imread('Clay_bit1.tiff')); %Name of binary image
BWoutline=bwperim(bw); %Extracts perimeter of the binary image

%% FRACTAL DIMENSIONS

%Fractal dimensions of the complete image
df_2D=hausDim(bw);% 2D Fractal dimension of the whole image using box counting 
dfhaus=dfhaus1(BWoutline);%Fractal dimension of the whole image based on the perimeter (Hausdorff dimension)
df_circle=dfcircle((bw));% 2D Fractal dimension of the whole image using concentric circles 
df_sandbox=dfsand(bw);% Fractal dimension found with sand boxing

% Fractal dimension based on the elipses that can be fitted to each blob in the image 
label=bwlabel((bw),8); 
blobmeasurements=regionprops(label,bw,'all');
numbofblobs=size(blobmeasurements,1);

phi=linspace(0,2*pi,50);
cosphi=cos(phi);
sinphi=sin(phi);

for k = 1 : numbofblobs           % Loop through all blobs.
	
	thisBlobsPixels = blobmeasurements(k).PixelIdxList;  % Get list of pixels in current blob.
	blobArea = blobmeasurements(k).Area;		% Get area.
	blobPerimeter = blobmeasurements(k).Perimeter;		% Get perimeter.
	blobCentroid = blobmeasurements(k).Centroid;		% Get centroid one at a time
	blobECD(k) = sqrt(4 * blobArea / pi);				% Compute ECD - Equivalent Circular Diameter.
    
    thisBlobsBoundingBox = blobmeasurements(k).BoundingBox;  % Get coordinates of the bounding box of the blob
    xbar = blobmeasurements(k).Centroid(1);
    ybar = blobmeasurements(k).Centroid(2);

    a = blobmeasurements(k).MajorAxisLength/2;
    b = blobmeasurements(k).MinorAxisLength/2;

    theta = pi*blobmeasurements(k).Orientation/180;
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;
    
    areal(k)=blobArea;
    majoral(k)=a;
end

logmajoral=log(majoral);
logareal=log(areal);

ft = fittype( 'poly1' );
[fitresult, gof] = fit( logmajoral', logareal', ft );
c11=fitresult.p1;

df_ellipse=c11(1);% Fractal dimension based on ellipses fitted to the blobs
rsqrdf1=gof.rsquare;


% Fractal dimensions of the infividual blobs

for k = 1 : numbofblobs 
    
    Pixels=blobmeasurements(k).BoundingBox;
    Pixels=[Pixels(1)-10 Pixels(2)-10 Pixels(3)+10 Pixels(4)+10];
    thisblob=ismember(label,k);
    J=(imcrop(thisblob,Pixels));
    J = bwareaopen(J,200,8);
    
    df2_individual(k)=hausDim(J);
    Joutline=bwperim(J);
    dfhaus_perim_ind(k)=dfhaus1(Joutline);
    df_circle_ind(k)=dfcircle((J));
    df_sand_ind(k)=dfsand(J);
    
end

% Average and standard deviation of the individual blobs

      df2_individual_1=nanmean(df2_individual);
      df2_individual_1_std=nanstd(df2_individual);
      dfhaus_perim_ind(find(dfhaus_perim_ind==Inf))=NaN;
      dfhaus_perim_ind(find(dfhaus_perim_ind==-Inf))=NaN;
      dfhaus_perim_ind_1=nanmean(dfhaus_perim_ind);
      dfhaus_perim_ind_1_std=nanstd(dfhaus_perim_ind);
      df_circle_ind_1=mean(nonzeros(df_circle_ind));
      df_circle_ind_1_std=std(nonzeros(df_circle_ind));
      df_sand_ind_1=mean(nonzeros(df_sand_ind));
      df_sand_ind_1_std=std(nonzeros(df_sand_ind));
