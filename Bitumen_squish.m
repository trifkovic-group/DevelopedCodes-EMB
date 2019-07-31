% 
function[Out1]=Bitumen_squish(z,bw_bitumen)
% Code written by Edna MB and Allen Shi
% Code for squishing the bitumen 3-D reconstruction by eliminating
% repetitive layers. Corrects Z position
% Last modification date: July 17th-2019
% z=100;

Raw=bw_bitumen;% Array defined in line 81 of Thresholding_Clay_Bitumen.m   
Raw = bwareaopen(Raw,200,26);% removes very small elements that probably come from noise, 200 means that elements with sizes lower than 200 are erased and 26 is the connectivity level.
CC = bwconncomp(Raw);% Finds connected components in the binary 3-D image.
RP = regionprops(CC,'All');% Measure the properties of the connected components
Out1 = zeros(size(Raw));%Predefines the output array with the squished elements
[X,Y,Z]=size(Raw);
sim=0.4; %similarity

for i=1:length(RP)%Evaluates each connected component found in line 10
 
RPC = abs(RP(i).BoundingBox(end,3)+RP(i).BoundingBox(end,6))/2; %finds the center of each component
%Auxiliar variables
count=0;
initial=1;
final=0;
A1=0;
B=0;
out_element=zeros(X,Y,Z);
P_element=[];
p_element=[];

if RP(i).PixelList(1,3)==1 && RP(i).PixelList(end,3)==z %if the object touches the top layer or bottom layer, don't move it
    Out1(RP(i).PixelIdxList)=1;

else
[x1pl,y1pl]=size(RP(i).PixelList); %gives you the size of the vector that stores the coordinates of each pixel that is part of each drop of bitumen

%% this section defines the coordinates of the pixels of each layer (Z)
for j=1:x1pl-1
A1=RP(i).PixelList(j,:); 
B=RP(i).PixelList(j+1,:); 
if A1(3)==B(3) 
else
count=count+1;
initial(count+1)=j+1;
final(count)=(j);
end
end

if size(initial)==1
final=[j+1];
else 
final=[final j+1];
end

Mat=[];
New_image=[];
comp=[];
Mat=[initial' final'];
New_image=RP(i).PixelList(1:final(1),:);
Out = zeros(size(Raw));
comp=RP(i).PixelList(initial(1):final(1),1:3);
count1= RP(i).PixelList(1,3);

mati=size(Mat);

%% This section eliminates the similar layers (They have to be different in at least dif% to be considered different)


for k=2:mati(1)
    
   A2=[];
   A2=RP(i).PixelList(initial(k):final(k),1:3);
   [X1A2,Y2A2]= size(A2);
   [X1comp,Y2comp]= size(comp);
   maxLength=max([X1A2,X1comp]); 
   A2(X1A2+1:maxLength,:) = 0;
   comp(X1comp+1:maxLength,:) = 0;
  [X3,Y3]=size(A2);

   if sum(sum(A2(:,1:2)==comp(:,1:2)))/(X3*(Y3-1))>sim
       comp=comp(any(comp,2),:);       
   else
    count1=count1+1;
       new_A1 = A2(any(A2,2),:);
       [X1,Y1]=size(new_A1);
       vec=ones(X1,1)*count1;
       new_A11=[new_A1(:,1:2) vec];
       New_image=[New_image;new_A11];
       comp=new_A1(any(new_A1,2),:); 
   end
end 

RPPI=0;

% The position of the reconstructed drop of bitumen only considers the top layer
RPPI = sub2ind(size(Out),New_image(:,1),New_image(:,2),New_image(:,3));

Out(RPPI) = 1;
cc2=regionprops(Out,'BoundingBox');
centroid=abs(cc2.BoundingBox(3)+cc2.BoundingBox(6))/2;

% Moves the centroid of the reconstructed drop to that of the original drop of bitumen

Size_New_image=size(New_image);

if New_image(1,3)==1  
RC_Plus=0;
elseif RP(i).PixelList(end,3)==z
        RC_move = (z-New_image(end,3));
        RC_Plus=[zeros(Size_New_image(1),2) ones(Size_New_image(1),1)*RC_move];
elseif New_image(1,3)~=1 && RP(i).PixelList(end,3)~=z
    RC_move =(RPC-centroid);
    RC_Plus=[zeros(Size_New_image(1),2) ones(Size_New_image(1),1)*RC_move];
end


New_image1=New_image+RC_Plus;


RPPI1 = sub2ind(size(Out1),New_image1(:,1),New_image1(:,2),round(New_image1(:,3)));
out_element(RPPI1)=1;
out_element=rot90(fliplr(out_element),1);

CC_element = bwconncomp(out_element);
P_element= regionprops(CC_element,'PixelIdxList');
p_element = P_element.PixelIdxList;
Out1(p_element) = 1;

end

end
end


