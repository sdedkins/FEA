function [Cor, Pos,sample,total,length_in_styecast] = Preprocessor( aspectab,aspectac,meshx_sample,meshx_sty,meshx,meshy,meshz,bcmode,styecastratio,nodes,nodex,nodey,nodez )
%Preprocessor - A function that returns matrices containing the coordinates
%and topology of the nodes. 
%   Detailed explanation goes here



% Define dimensions in arbitrary units normalising the length along the x
% direction to 100

samplea=100;
sampleb=samplea/aspectab;
samplec=samplea/aspectac;

length_in_styecast=styecastratio*min(sampleb,samplec);



if bcmode==1
    totala=samplea;
    totalb=sampleb;
    totalc=samplec;
    
elseif bcmode==2||bcmode==3
    totala=samplea+2*length_in_styecast;
    totalb=sampleb;
    totalc=samplec;
else error('Boundary Condition Mode Not Recognised')
    
end

% Put sample and total dimensions into vectors for returning to main
% program

sample=[samplea,sampleb,samplec];
total=[totala,totalb,totalc];



%Mesh error control
if meshx <1 || meshy <1 || meshz <1;
    error('Open Mesh - Redefine number of elements')
else if meshx >100 || meshy>15|| meshz >15;   
    display('Mesh is Large - This may take some time')
    end
end


% Generate positions of nodes. Nodes are equally spaced along each direction (but with
% different spacings for different directions). For bcmode = 2 or 3 the
% extra length of material embedded in stycast is included


if bcmode==1
    
xcor=linspace(0,total(1),nodex);
ycor=linspace(0,total(2),nodey);
zcor=linspace(0,total(3),nodez);

elseif bcmode==2 || bcmode==3;
    
    xcor=[linspace(0,length_in_styecast,meshx_sty+1),...
        linspace(length_in_styecast,length_in_styecast+sample(1),meshx_sample+1),...
        linspace(length_in_styecast+sample(1),total(1),meshx_sty+1)];
    xcor=unique(xcor);
    
    ycor=linspace(0,total(2),nodey);
zcor=linspace(0,total(3),nodez);

else
    
    error('Unrecognised Boundary Condition Mode');
end
    
% Generate array Cor which holds are list of the node co-ordinates for
% each node number i.e Cor(n,:)=[xcor for node n, ycor for node n, zcor for
% node n]
            
l=1;
for k=1:nodez
    for j=1:nodey
        for i=1:nodex
       
            Cor(l,:)=[xcor(i),ycor(j),zcor(k)];
            l=l+1;
        end
    end
end

clear i j k l ;


% So far the mesh has been generated for x,y,z all >= 0. If ends are
% embedded in stycast this bit of code displaces the whole mesh in the
% negative x direction so that the sample outside of the stycast starts at
% x=0 and ends at x=100.

if bcmode==2 || bcmode==3
    
   Cor(:,1)=Cor(:,1)-length_in_styecast;
end  



%Calculate the Matrix Pos which is a list of which nodes belond to which elements      
 %Pos( Element No,:)=[Shell Element Node Number]         
no=1;             
for k=1:meshz
    for i=1:meshy;
        for j=1:meshx;
            if  i < meshy+1 && j<meshx+1 ;
                
 Pos(no,:)=[(meshx+1)*(i-1)+j+(k-1)*(meshx+1)*(meshy+1)   ...
            (meshx+1)*(i-1)+j+1+(k-1)*(meshx+1)*(meshy+1) ... 
            (meshx+1)*(i)+j+1+(k-1)*(meshx+1)*(meshy+1)   ... 
            (meshx+1)*(i)+j+(k-1)*(meshx+1)*(meshy+1)     ... 
            (meshx+1)*(i-1)+j+(k)*(meshx+1)*(meshy+1)     ...
            (meshx+1)*(i-1)+j+1+(k)*(meshx+1)*(meshy+1)   ...
            (meshx+1)*(i)+j+1+(k)*(meshx+1)*(meshy+1)     ...
            (meshx+1)*(i)+j+(k)*(meshx+1)*(meshy+1)];

            end
            no=no+1;    
        end
    end
end%________________________________________|

clear i j k no;






end

