function [ SX,SY,SZ,SSTRAIN ] = Sliceform( dimensions,strain,xcen,ycen,zcen )
%Turns data outputed from FEA into form required by Slice function
%   Detailed explanation goes here
xelements=dimensions(4);
yelements=dimensions(5);
zelements=dimensions(6);



%%____Define Volume Matrices for Slice Function______%%

SX=zeros(yelements,xelements,zelements);


for i=1:zelements
    for j=1:yelements
        

SX(j,:,i)=xcen((1+xelements*(j-1)+xelements*yelements*(i-1)):xelements*(j-1)+xelements*yelements*(i-1)+xelements)';

    end
    
end

clear i j

for i=1:zelements
    for j=1:yelements
        

SY(j,:,i)=ycen((1+xelements*(j-1)+xelements*yelements*(i-1)):xelements*(j-1)+xelements*yelements*(i-1)+xelements)';

    end
    
end

for i=1:zelements
    for j=1:yelements
        

SZ(j,:,i)=zcen((1+xelements*(j-1)+xelements*yelements*(i-1)):xelements*(j-1)+xelements*yelements*(i-1)+xelements)';

    end
    
end

for k=1:6
    for i=1:zelements
     for j=1:yelements
        

SSTRAIN(j,:,i,k)=strain(k,(1+xelements*(j-1)+xelements*yelements*(i-1)):(xelements*(j-1)+xelements*yelements*(i-1)+xelements))';

     end
    
    end
end



end

