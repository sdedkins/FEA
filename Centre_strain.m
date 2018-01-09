function [ CENTRE_STRAIN ] = Centre_strain( H_DATA,aspect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% Determine co-ordinate of centre point of sample for each simulation

CENTRE_STRAIN = struct('x',{},'y',{},'z',{},'exx',{},'eyy',{},'ezz',{}); % Define empty structure array
%to hold centre co-ordinates

for i=1:length(aspect)
    for bcmode=1:length(H_DATA(1,:,1))
        for matched=1:length(H_DATA(1,1,:))
            
            
            SX=H_DATA(i, bcmode,matched).SX;
            SY=H_DATA(i, bcmode,matched).SY;
            SZ=H_DATA(i, bcmode,matched).SZ;
            SSTRAIN=H_DATA(i, bcmode,matched).SSTRAIN;
            xi=H_DATA(i,bcmode,matched).dimensions(1)/2;
            yi= H_DATA(i,bcmode,matched).dimensions(2)/2;
            zi= H_DATA(i,bcmode,matched).dimensions(3)/2;
            
            
            CENTRE_STRAIN(i,bcmode,matched).x=xi;
            
            CENTRE_STRAIN(i,bcmode,matched).y= yi;
            
            CENTRE_STRAIN(i,bcmode,matched).z=zi;
            
            CENTRE_STRAIN(i,bcmode,matched).exx=interp3(SX,SY,SZ,SSTRAIN(:,:,:,1),xi,yi,zi);
            
            CENTRE_STRAIN(i,bcmode,matched).eyy=interp3(SX,SY,SZ,SSTRAIN(:,:,:,2),xi,yi,zi);
            
            CENTRE_STRAIN(i,bcmode,matched).ezz=interp3(SX,SY,SZ,SSTRAIN(:,:,:,3),xi,yi,zi);
            
        end
    end
end








end

