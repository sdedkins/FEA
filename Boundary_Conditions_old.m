function [ Re,R,dofConstraints,P ] = Boundary_Conditions( Pos,ex, ey, ez,nodes,nodex,nodey,nodez,meshx,meshy,meshz,nel,ndof,bcmode,sample,total)
%Boundary_Conditions Summary of this function goes here
%   Detailed explanation goes here





%_____THIS PART OF CODE ASSIGNS dof#s to nodes____%


%Re(Node number,:)=[xdof# ydof# zdof#]

Re=zeros(nodes,3);
k=0;
for i=1:nodes;
    for j=1:3;
       
            k = k +1;
            Re(i,j) =k ;
        
    end
end

clear i j k



%Calculate R: A Matrix that tells you for a given element number what the dof#s for its nodes are 

% R(element number,direction)= degrees of freedom of first node in element, degrees of
%freedom of node 2nd node in element,...
for i=1:nel
R(i,:)=[Re(Pos(i,1),:) Re(Pos(i,2),:) Re(Pos(i,3),:) Re(Pos(i,4),:)...
        Re(Pos(i,5),:) Re(Pos(i,6),:) Re(Pos(i,7),:) Re(Pos(i,8),:)];
end



%____THIS PART OF CODE ASSIGNS CONSTRAINTS TO THE APPROPRIATE dofs________%

disp('Calculating Constraints')


% Calculate displacements of boundary nodes that correspond to the enforced
% strain

ux=sample(1)*ex/2;
uy=sample(2)*ey/2;
uz=sample(3)*ez/2;

if bcmode==1
    

        Constraintsx=zeros(2*nodey*nodez,2);


        % Apply x constraint to the end faces
        for j=1:nodez
            for i=1:nodey

                Constraintsx(i+(j-1)*nodey,:)=[1+(i-1)*nodex+(j-1)*nodex*nodey,ux];
            end
        end 
        clear i j



        % now apply y and z constraints to the edges of the end faces

        Constraintsz=zeros(4*nodey,2);
        Constraintsy=zeros(4*nodez,2);

        Constraintsz(1:nodey,1)=1:nodex:nodey*nodex;
        Constraintsz(1:nodey,2)=uz;
        Constraintsz(nodey+1:(2*nodey),:)=Constraintsz(1:nodey,:)+repmat([(meshz)*(nodex)*(nodey),-2*uz],nodey,1);
        Constraintsy(1:(nodez),:)=[(1:(nodex)*(nodey):(nodex)*(nodez)*(nodey))', repmat(uy,nodez,1)];
        Constraintsy((nodez+1):(2*(nodez)),:)= Constraintsy(1:(nodez),:) + repmat([(meshy)*(nodex), -2*uy],nodez,1);


        %now copy to otherside of block

        Constraintsx((nodez)*(nodey)+1:end,:)=Constraintsx(1:(nodez)*(nodey),:)+repmat([(meshx),-2*ux],(nodez)*(nodey),1);
        Constraintsz(2*(nodey)+1:end,:)=Constraintsz(1:(2*(nodey)),:)+repmat([(meshx),0],2*(nodey),1);
        Constraintsy(2*(nodez)+1:end,:)=Constraintsy(1:(2*(nodez)),:)+repmat([(meshx),0],2*(nodez),1);


        %Copy into matrix that lists [dof# constrained value]

        nocx=size(Constraintsx,1);
        nocy=size(Constraintsy,1);
        nocz=size(Constraintsz,1);

        dofConstraints(1:nocx,:)=[Re(Constraintsx(:,1),1), Constraintsx(:,2)];
        dofConstraints((nocx+1):(nocx+nocy),:)=[Re(Constraintsy(:,1),2), Constraintsy(:,2)];
        dofConstraints((nocx+nocy+1):(nocx+nocy+nocz),:)=[Re(Constraintsz(:,1),3), Constraintsz(:,2)];

elseif bcmode==2
        
    %Apply x constraints to end faces
    
    Constraintsx=(Pos(1:nelx:end,
        
        
        
else
        error('bcmode has taken illegal value')
        
end


%%______________________________________THIS PART OF CODE DEFINES THE
%%LOADS ON THE SYSTEM___________________________%%

P=zeros(1,ndof);



end

