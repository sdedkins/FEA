function [ Re,R,dofConstraints_unique,P ] = ...
Boundary_Conditions( Pos,Cor,ex, ey, ez,nodes,nel,ndof,bcmode,sample,total)

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



%Calculate R: A Matrix that tells you for a given element number 
% what the dof#s for its nodes are 

% R(element number,direction)= [degrees of freedom of first node in element
%,degrees of freedom of node 2nd node in element,...]
for i=1:nel
R(i,:)=[Re(Pos(i,1),:) Re(Pos(i,2),:) Re(Pos(i,3),:) Re(Pos(i,4),:)...
        Re(Pos(i,5),:) Re(Pos(i,6),:) Re(Pos(i,7),:) Re(Pos(i,8),:)];
end



%____THIS PART OF CODE ASSIGNS CONSTRAINTS TO THE APPROPRIATE dofs________%




% Calculate displacements of boundary nodes that correspond to the enforced
% strain

ux=sample(1)*ex/2;
uy=sample(2)*ey/2;
uz=sample(3)*ez/2;



if bcmode==1||bcmode==2
    
    

        % find degrees of freedom corresponding to the x co-ordinate of 
        %the nodes on the end faces
        constrained_end1= find(Cor(:,1)==min(Cor(:,1)));
        constrained_end2= find(Cor(:,1)==max(Cor(:,1))) ;
        % Apply the x displacement constraint to these
        dofConstraints_end=...
            [Re(constrained_end1,1),repmat(-ux,length(constrained_end1),1);
            Re(constrained_end2,1),repmat(+ux,length(constrained_end2),1)];



        % Now find all the dof#s corresponding the the y coordinate of nodes
        %on the faces perpendicular to the y-direction that are constrained
        %by stycast
        constrainedy1=find(Cor(:,2)==0 & (Cor(:,1)<=0 | Cor(:,1)>=sample(1)));
        constrainedy2=find(Cor(:,2)==total(2) & (Cor(:,1)<=0 | Cor(:,1)>=sample(1)));
        %Apply the y-displacement constraint to these
        dofConstraintsy=[Re(constrainedy1,2),repmat(-uy,length(constrainedy1),1);
                         Re(constrainedy2,2),repmat(+uy,length(constrainedy2),1)];

         % Now find all the dof#s corresponding the the z coordinate of nodes
        %on the faces perpendicular to the z-direction that are constrained
        %by stycast             
        constrainedz1=find(Cor(:,3)==0 & (Cor(:,1)<=0 | Cor(:,1)>=sample(1)));
        constrainedz2=find(Cor(:,3)==total(3) & (Cor(:,1)<=0 | Cor(:,1)>=sample(1)));
        %Apply the z-displacement constraint to these
        dofConstraintsz=[Re(constrainedz1,3),repmat(-uz,length(constrainedz1),1);
                         Re(constrainedz2,3),repmat(+uz,length(constrainedz2),1)];
        %

        
         % Now find all the dof#s corresponding the the x coordinate of nodes
        %on the outside surfaces that are not the end faces       
        constrainedx1=find((Cor(:,1)<=0)&(Cor(:,2)==0|Cor(:,2)==total(2)...
            |Cor(:,3)==0|Cor(:,3)==total(3)));

        constrainedx2=find((Cor(:,1)>=sample(1))&(Cor(:,2)==0|Cor(:,2)==total(2)...
            |Cor(:,3)==0|Cor(:,3)==total(3)));
        %Apply x-displacment condition
        dofConstraintsx=[Re(constrainedx1,1),repmat(-ux,length(constrainedx1),1);
                          Re(constrainedx2,1),repmat(+ux,length(constrainedx2),1)];
        % Concantenate

        dofConstraints=[dofConstraints_end;dofConstraintsy;dofConstraintsz;dofConstraintsx];


        % Make unique (get rid of duplicate entries).

        dofConstraints_unique=unique(dofConstraints, 'rows');
        
        
        %Phew...glad thats over. Now time to do the same all over again
        % but for the case where only the bottom face in costrained by
        % stycast

elseif bcmode==3
    
       % Find all dof# corresponding to y cor of nodes on bottom of bit in
       % stycast and apply displacement to them
        constrainedy1=find(Cor(:,2)==0 & (Cor(:,1)<=0 | Cor(:,1)>=sample(1))&(Cor(:,3)==0));
        constrainedy2=find(Cor(:,2)==total(2) & (Cor(:,1)<=0 | Cor(:,1)>=sample(1))&(Cor(:,3)==0));
        dofConstraintsy=[Re(constrainedy1,2),repmat(-uy,length(constrainedy1),1);
                         Re(constrainedy2,2),repmat(+uy,length(constrainedy2),1)];
                     
        % Find all dof# corresponding to x cor of nodes on bottom of bit in
       % stycast and apply displacement to them
        constrainedx1=find((Cor(:,1)<=0)&Cor(:,3)==0);

        constrainedx2=find((Cor(:,1)>=sample(1))&Cor(:,3)==0);

        dofConstraintsx=[Re(constrainedx1,1),repmat(-ux,length(constrainedx1),1);
                          Re(constrainedx2,1),repmat(+ux,length(constrainedx2),1)];    
         
                      
        % Constrain bottom face of bit in stycast to have 0 displacement 
        %in z direction. Required if not to be underconstrained.                
        constrainedz1=find(Cor(:,1)<=0&Cor(:,3)==0);
        constrainedz2=find(Cor(:,1)>=sample(1)&Cor(:,3)==0);
        dofConstraintsz=[Re(constrainedz1,1),repmat(0,length(constrainedz1),1);
                          Re(constrainedz2,1),repmat(0,length(constrainedz2),1)];
        
        % concantenate

        dofConstraints=[dofConstraintsy;dofConstraintsx;dofConstraintsz]; 
        
          % make unique

        dofConstraints_unique=unique(dofConstraints, 'rows');
        
else
    error('Unrecognised Boundary Condition Value');
    
end
        
        


%%______________________________________THIS PART OF CODE DEFINES THE
%%LOADS ON THE SYSTEM___________________________%%

P=zeros(1,ndof);



end

