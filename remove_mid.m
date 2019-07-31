function [Out1] = remove_mid(AreaOpenedIm_i1)
%remove mid_removes objects that are not touching top or bottom.  

Raw=AreaOpenedIm_i1;% Final thresholded images
Raw = bwareaopen(Raw,3000,26);% Removes very small elements that probably come from noise, 200 means that elements with sizes lower than 200 are erased and 26 is the connectivity level.
CC = bwconncomp(Raw);% Finds connected components in the binary 3-D image.
RP = regionprops(CC,'All');% Measure the properties of the connected components
Out1 = zeros(size(Raw));%Predefines the output array with the squished elements
[X,Y,Z]=size(Raw);

for i=1:length(RP)%Evaluates each connected component found in line 6
if RP(i).PixelList(1,3)==1 || RP(i).PixelList(end,3)==Z
    Out1(RP(i).PixelIdxList)=1;    
elseif RP(i).PixelList(1,3)~=1 && RP(i).PixelList(end,3)~=Z

end
end

end
