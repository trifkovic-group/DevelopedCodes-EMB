
function df_circle = dfcircle(bw)
% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
[y x]=size(bw);

i=linspace(1,y/2,10);

for k=length(i):-1:1

    %%creates de circle
imageSizeX = x;
imageSizeY = y;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = x/2;
centerY = y/2;
radius =i(k);
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
% circlePixels is a 2D "logical" array.
% Now, display it.
% figure, image(circlePixels) 
colormap([0 0 0; 1 1 1]);
% title('Binary image of a circle');

bw1=circlePixels.*bw;
% figure, imshow(bw1)

[I J]=find(bw1==1);
%area(k)=length(I);
area(k)=bwarea(bw1);
radius_k(k)=radius;

end

zl=find(area);
log_area(zl)=log(area(zl));
log_radius(zl)=log(radius_k(zl));

if length(log_area)<2
    df_circle=0;
    
else


ft = fittype( 'poly1' );
[fitresult, gof] = fit(log_radius',log_area', ft );
c12=fitresult.p1;

df_circle=c12(1);
rsqrdf2=gof.rsquare;

   
    
end



end