function df_sand2 = dfsand(bw)
% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
[y x]=size(bw);

i=linspace(1,log(y/2),10);

for k=length(i):-1:1

imageSizeX = x;
imageSizeY = y;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = x/2;
centerY = y/2;
width =exp(i(k))*2;
height=exp(i(k))*2;
xleft=centerX-width/2;
ybottom=centerY-height/2;
squarePixels =[xleft, ybottom, width, height];

% circlePixels is a 2D "logical" array.
% Now, display it.
% figure, image(circlePixels) 

% title('Binary image of a circle');

bw1=imcrop(bw,squarePixels);
% figure, imshow(bw1)

[I J]=find(bw1==1);
area(k)=length(I);
%area(k)=bwarea(bw1);
radius_k(k)=width;

end

zl=find(area);
log_area(zl)=log(area(zl));
%zl=find(radius_k);
log_radius(zl)=log(radius_k(zl));

if length(log_area)<2
    df_sand2=0;
    
else


ft = fittype( 'poly1' );
[fitresult, gof] = fit(log_radius',log_area', ft );
c12=fitresult.p1;

df_sand2=c12(1);
rsqrdf2=gof.rsquare;

% plot(log_radius(2:end)',log_area(2:end)');
    
    
end



end