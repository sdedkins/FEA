function [ D , sparseKsis] = Solver( Cor,Pos,R,E,dofConstraints,P,nel,ndof)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%_____________________________GAUSS LEGENDRE CONSTANTS
%format long;
l=1;
gpoint(1)=-0.577350269189626;
gpoint(2)=0.577350269189626;
for k=1:2
    for j=1:2
        for i=1:2
W(l,:)=[gpoint(i),gpoint(j),gpoint(k)];
l=l+1;
        end
    end    
end
%_________________________________________________|
                  
clear l i j k;



tic
% Preallocate Variables

pxn=zeros(8,2);
pyn=zeros(8,2);
pzn=zeros(8,2);
A=zeros(6,24,nel);
K(24,24,nel)=0;


%_____________________________LOCAL AXIS ELEMENT STIFFNESS MATRIX
for s=1:nel;%total element number


    
%Moving global cartesian nodes to element nodes
for u=1:8
pxn(u,:) =Cor(Pos(s,u),1);  %x(u)
pyn(u,:) =Cor(Pos(s,u),2);  %y(u)
pzn(u,:) =Cor(Pos(s,u),3);  %z(u)
end



for i=1:8;
        %[i j k s] Numerical integration numerator
        %Parametric parameters is connection homogen parameters

e=W(i,1);
n=W(i,2);
J=W(i,3);


%____FUNCTION____:Element shape functions partial derivatives.
[Jacobi,Hx,Hy,Hz]=AconnectH8(pxn,pyn,pzn,e,n,J);

%[A]=[Equbilibrium differantial matrix].[shell shape matrix]
%[A]=[d].[N]

for f=1:8;
A(:,3*(f-1)+1,s)= [Hx(f);   0   ;   0  ;  Hy(f) ;   0   ;  Hz(f)]; 
A(:,3*(f-1)+2,s)= [ 0   ; Hy(f) ;   0  ;  Hx(f) ;  Hz(f);   0   ]; 
A(:,3*(f-1)+3,s)= [ 0   ;   0   ; Hz(f);   0    ;  Hy(f);  Hx(f)]; 
end

              K(:,:,s)=K(:,:,s)+det(Jacobi)*A(:,:,s)'*E*A(:,:,s);

%fprintf ('Element integrate number = % .2f \n',s);
     
end

end%____________________________|



disp('Calculating System Stiffness Matrix')


Ksis=zeros(ndof,ndof);
%___________________________SYSTEM STIFFNESS MATRIX
for n=1:nel;
        for sat=1:24;
            for sut=1:24;
               
                    Ksis(R(n,sut),R(n,sat))=Ksis(R(n,sut),R(n,sat)) + K(sat,sut,n);
              
            end%
        end%
end%_______________________________________________|
%Ksis:Global system stiffness matrix
clear n sat sut 




% Apply Essential Boundary Conditions

bcwt=trace(Ksis)/24;

P =P'- Ksis(:,dofConstraints(:,1))*dofConstraints(:,2);
Ksis(:,dofConstraints(:,1)) = 0;
Ksis(dofConstraints(:,1),:) = 0;
Ksis(dofConstraints(:,1),dofConstraints(:,1)) = bcwt*speye(length(dofConstraints(:,1)));
P(dofConstraints(:,1)) = bcwt*dofConstraints(:,2);


sparseKsis=sparse(Ksis);

clear Ksis;
%Ksis:System global stiffness matrix
%______________________________SYSTEM DISPLACEMENT
%equation=size(Ksis);
%if equation(1)~=rank(Ksis)
 %   display('This system stiffness matrix is badly scaled')
%    R
  %  error('Control system support boundary conditions')
%else
   
    D = sparseKsis\P;

clear equation Ku





end

