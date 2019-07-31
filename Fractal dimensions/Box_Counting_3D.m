function [Df_3d] = Box_Counting_3D(resxy,resz,Z,Clay_bit)
%Code written by Edna MB
%This program determines the 3D-fractal dimension of a structure
%Last modification date: February 20th-2019

P=resz*(1:1:Z);
R=resxy*(1:1:Z);

% If the resolution xy and z is different, this section equalizes resolutions
k=0;
for i=1:length(R)
   for j=1:length(P)-1 
   if R(i)>P(j) && R(i)<P(j+1)
       k=k+1;
       BW_nuevo(:,:,i)=Clay_bit(:,:,j);
       BW_nuevog(:,:,i)=Clay_bit(:,:,j);
       tra(k)=j;
   end
      end
end

e2=(BW_nuevo);
Nx=size(e2,1);
Ny=size(e2,2);
Nz=size(e2,3);

al=1:1:Nz;

%Box-counting 
for np = 1:length(al);
  
    ali=(al(np));
    numBlocks=floor(Nz/ali);
    
    sizeBlocks_x = floor(Nz./numBlocks);
    sizeBlocks_y = floor(Nz./numBlocks);
    sizeBlocks_z= floor(Nz./numBlocks);

    numBlocksx=floor(Nx/sizeBlocks_x);
    numBlocksy=floor(Ny/sizeBlocks_y);
    numBlocksz=floor(Nz/sizeBlocks_z);
    
    flag = zeros(numBlocksx,numBlocksy,numBlocksz);
    for l = 1:numBlocksx
        for j = 1:numBlocksy
            for k=1:numBlocksz
            xStart = (l-1)*sizeBlocks_x +1;
            xEnd   = l*sizeBlocks_x;
        
            yStart = (j-1)*sizeBlocks_y + 1;
            yEnd   = j*sizeBlocks_y;
            
            zStart = (k-1)*sizeBlocks_z + 1;
            zEnd   = k*sizeBlocks_z;
            
            block = e2(xStart:xEnd, yStart:yEnd, zStart:zEnd);
        
            flag(l,j,k) = any(block(:)); %mark this if ANY part of block is true
            end
        end
    end
    boxCount = nnz(flag);
    resolution(np)=numBlocksz;
    %resolution(numBlocks)=sizeBlocks_x;
    table(np)    = boxCount;
end

x1=log((resolution));
y1=log((table));
x1=[x1(1:7) x1(10) x1(14)];
y1=[y1(1:7) y1(10) y1(14)];


p2       = polyfit(x1,y1,1);
BestFit2 = polyval(p2,x1);
Df_3d = p2(:,1);

end

