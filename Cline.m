function [ clinex, clines ] = Cline( SX,SY,SZ,SSTRAIN )
% A function that returns data along the centre line of the volume for
% plotting
%   Detailed explanation goes here

xelements=size(SX,2);
yelements=size(SX,1);
zelements=size(SY,3);

if mod(yelements,2)==1
    
    a(:,:)=SX(0.5*(yelements+1),:,:);
astrain(:,:,:)=SSTRAIN(0.5*(yelements+1),:,:,:);
astrain=squeeze(astrain);
    else
    a(:,:)=0.5*(SX(0.5*yelements,:,:)+SX(0.5*yelements+1,:,:));
    astrain=0.5*(SSTRAIN(0.5*yelements,:,:,:)+SSTRAIN(0.5*yelements+1,:,:,:));
    astrain=squeeze(astrain);
end

if mod(zelements,2)==1
    
    clinex(:)=a(:,0.5*(zelements+1));
   clines(:,:)=astrain(:,0.5*(zelements+1),:);
    clines=squeeze(clines);
else
    
    clinex=0.5*(a(:,0.5*zelements)+a(:,0.5*zelements+1));
    clines(:,:)=0.5*(astrain(:,0.5*zelements,:)+astrain(:,0.5*zelements+1,:));
    clines=squeeze(clines);
end
end

