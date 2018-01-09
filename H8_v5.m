% _________________________________________________________________________
%| 8NODE-24DOF HEXAHEDRAL SOLID ISOPARAMETRIC   F.E.M.     [A.Ö] (C)(R)    |
%|##########################################                               |
%| MECHANISM           | Edge       Upperside  4X(1Node-3Dof)[ux,vy,wz]    |
%|                     |            Underside  4X(1Node-3Dof)[ux,vy,wz]    |
%|                     | Middle     node      N/A                          |
%|                     | Centroidal node      N/A                          |
%|---------------------|---------------------------------------------------|
%| ANALYSIS            | Shape function:Parametric (Pascal triangle)       |
%|                     | Stifness      :Topology+Accumulate method         |
%|                     | Solve Equation:Cholesky-[L][D][u]                 |
%|-------------------------------------------------------------------------|
%|SUB FUNCTIONS                                                            |
%|1- AconnectH8(pxn,pyn,pzn,e,n,J)                    [AconnectH8.m]       |
%|2- threedimensionalstress(sigmax,sigmay,sigmaz,tauxy,tauyz,tauxz,1)      |
%|-------------------------------------------------------------------------|
%| Maxima 5.9.0.9beta2 http://maxima.sourceforge.net                       |
%| Using Lisp Kyoto Common Lisp GCL 2.6.3 (aka GCL)                        |
%| Distributed under the GNU Public License. See the file COPYING.         |
%| Dedicated to the memory of William Schelter.                            |
%| This is a development version of Maxima. The function bug_report()      |
%| provides bug reporting information.                                     |
%|_________________________________________________________________________|

%Shape function terms
%                         1
%                      x  y  z
%                  x*y   y*z  x*z
%                        y*z
%

%"ELEMENT PARAMETRIC SHAPE FUNCTIONS"
%N1:1/8*(1-2*x/a)*(1-2*y/b)*(1-2*z/c)
%N2:1/8*(1+2*x/a)*(1-2*y/b)*(1-2*z/c)
%N3:1/8*(1+2*x/a)*(1+2*y/b)*(1-2*z/c)
%N4:1/8*(1-2*x/a)*(1+2*y/b)*(1-2*z/c)
%N5:1/8*(1-2*x/a)*(1-2*y/b)*(1+2*z/c)
%N6:1/8*(1+2*x/a)*(1-2*y/b)*(1+2*z/c)
%N7:1/8*(1+2*x/a)*(1+2*y/b)*(1+2*z/c)
%N8:1/8*(1-2*x/a)*(1+2*y/b)*(1+2*z/c)

clc
clear 

tic
disp('Please wait Programme Running')



                                     %##################INPUT DATA###############
                                     
%####################### AUTOMATIC MESH DATA
totala  = 75;        %"x" axis direction total length  (m)
totalb  = 10;        %"y" axis direction total width   (m)
totalc  = 5;        %"z" axis direction total heigth  (m)
meshx   = 50;           %"x" axis direction mesh number
meshy   = 10;           %"y" axis direction mesh number
meshz   = 10;           %"z" axis direction segment value
segment = 2;           %Stress "z" axis directions values
ex=0;
ey=6.62e-4;
ez=1.83e-3;
E11=2.32;
E33=2.08;
E23=0.71;
E12=1.06;
E44=0.657;
E66=0.612;
%#######################

gridspacex = totala/meshx;  % along to "x" axis direction mesh space
gridspacey = totalb/meshy;  % along to "y" axis direction mesh space
gridspacez = totalc/meshz;  % along to "z" axis direction mesh space


disp('Defining Grid')

%Mesh error control
if meshx <1 || meshy <1 || meshz <1;
    display('Open dimension system mesh and this mesh value > 1')
%    [meshx meshy]
    error('mesh values')
else if meshx >10 || meshy>10 || meshz >10;   
    display('system analysis mesh is big')
    end
end

gridspacex1=0.75;
gridspacex2=4.5;


%_____________________________________GLOBAL COORDINATES
dim=1;
for ars=0:meshz
    for sut=0:meshy 
        x=-gridspacex1;
         for sat=0:meshx %Element edge node reference axis
          if sat <= 20
             x=x+gridspacex1;
          elseif sat >20 && sat <= 30
              x=x+gridspacex2;
          else
              x=x+gridspacex1;
          end    
          y=sut*gridspacey;
          z=ars*gridspacez;
%Global axis system nodes position
          Cor(dim,:)=[x y z];
         %Cor(dim,1)=x;    %"x" global coordinate value 
         %Cor(dim,2)=y;    %"y" global coordinate value 
         %Cor(dim,3)=z;
         dim=dim+1;    
        end
    end
end%___________________________________________|
clear x y z sat sut ars dim;
Node=size(Cor,1);                     %Elements Number

nodex=meshx+1;
nodey=meshy+1;
nodez=meshz+1;


%____________________________________POSITION MATRIX        
no=1;       %Pos(Shell Element No,:)=[Shell Element Node Number]        
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
No=size(Pos,1);   %Plane Elements Number

clear i j k;

nnodes=nodex*nodey*nodez;
ndof=3*nnodes;







%_____THIS PART OF CODE ASSIGNS dof#s to nodes____%


%Re(Node number,:)=[xdof# ydof# zdof#]


Re=zeros(Node,3);
k=0;
for i=1:Node;
    for j=1:3;
       
            k = k +1;
            Re(i,j) =k ;
        
    end
end


%_________________________________MODAL DISPLACEMENTS   % R(element number,
%direction)= degrees of freedom of first node in element, degrees of
%freedom of node 2nd node in element,...
for i=1:No
R(i,:)=[Re(Pos(i,1),:) Re(Pos(i,2),:) Re(Pos(i,3),:) Re(Pos(i,4),:)...
        Re(Pos(i,5),:) Re(Pos(i,6),:) Re(Pos(i,7),:) Re(Pos(i,8),:)];
end




%____THIS PART OF CODE ASSIGNS CONSTRAINTS TO THE APPROPRIATE dofs________%

disp('Calculating Constraints')


% Define system constraints

ux=totala*ex/2;
uy=totalb*ey/2;
uz=totalc*ez/2;

Constraintsx=zeros(2*(meshy+1)*(meshz+1),2);

for j=1:(meshz+1)
    for i=1:(meshy+1)

        Constraintsx(i+(j-1)*(meshy+1),:)=[1+(i-1)*(meshx+1)+(j-1)*(meshx+1)*(meshy+1),ux];
    end
end 
clear i j



% now apply y and z constraints to the edges

Constraintsz=zeros(4*(meshy+1),2);
Constraintsy=zeros(4*(meshz+1),2);

Constraintsz(1:(meshy+1),1)=[1:(meshx+1):((meshy+1)*(meshx+1))];
Constraintsz(1:(meshy+1),2)=uz;
Constraintsz(meshy+2:(2*(meshy+1)),:)=Constraintsz(1:(meshy+1),:)+repmat([(meshz)*(meshx+1)*(meshy+1),-2*uz],meshy+1,1);
Constraintsy(1:(meshz+1),:)=[(1:(meshx+1)*(meshy+1):(meshx+1)*(meshz+1)*(meshy+1))', repmat(uy,meshz+1,1)];
Constraintsy((meshz+2):(2*(meshz+1)),:)= Constraintsy(1:(meshz+1),:) + repmat([(meshy)*(meshx+1), -2*uy],meshz+1,1);


%now copy to otherside of block

Constraintsx((meshz+1)*(meshy+1)+1:end,:)=Constraintsx(1:(meshz+1)*(meshy+1),:)+repmat([(meshx),-2*ux],(meshz+1)*(meshy+1),1);
Constraintsz(2*(meshy+1)+1:end,:)=Constraintsz(1:(2*(meshy+1)),:)+repmat([(meshx),0],2*(meshy+1),1);
Constraintsy(2*(meshz+1)+1:end,:)=Constraintsy(1:(2*(meshz+1)),:)+repmat([(meshx),0],2*(meshz+1),1);


%Copy into matrix that lists [dof# constrained value]

nocx=size(Constraintsx,1);
nocy=size(Constraintsy,1);
nocz=size(Constraintsz,1);

dofConstraints(1:nocx,:)=[Re(Constraintsx(:,1),1), Constraintsx(:,2)];
dofConstraints((nocx+1):(nocx+nocy),:)=[Re(Constraintsy(:,1),2), Constraintsy(:,2)];
dofConstraints((nocx+nocy+1):(nocx+nocy+nocz),:)=[Re(Constraintsz(:,1),3), Constraintsz(:,2)];


%%______________________________________THIS PART OF CODE DEFINES THE
%%LOADS ON THE SYSTEM___________________________%%

P=zeros(1,ndof);





%"ELEMENT ELASTICITY MATRIX [C]";



C=zeros(6,6);

C(1,:)=[E11,E12,E23,0,0,0];
C(2,:)=[E12,E11,E23,0,0,0];
C(3,:)=[E23,E23,E33,0,0,0];
C(4,:)=[0,0,0,E44,0,0];
C(5,:)=[0,0,0,0,E44,0];
C(6,:)=[0,0,0,0,0,E66];

C=C.*1e11;



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


disp('Calculating Element Stiffness Matrix')
tic
% Preallocate Variables

pxn=zeros(8,2);
pyn=zeros(8,2);
pzn=zeros(8,2);
A=zeros(6,24,No);
K(24,24,No)=0;


%_____________________________LOCAL AXIS ELEMENT STIFFNESS MATRIX
for s=1:No;%total element number


    
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

              K(:,:,s)=K(:,:,s)+det(Jacobi)*A(:,:,s)'*C*A(:,:,s);

%fprintf ('Element integrate number = % .2f \n',s);
     
end

end%____________________________|



disp('Calculating System Stiffness Matrix')


Ksis=zeros(3*nnodes,3*nnodes);
%___________________________SYSTEM STIFFNESS MATRIX
for n=1:No;
        for sat=1:24;
            for sut=1:24;
               if (R(n,sat)~=0)
                  if (R(n,sut)~=0);
                    Ksis(R(n,sut),R(n,sat))=Ksis(R(n,sut),R(n,sat)) + K(sat,sut,n);
                  end%    
               end%
            end%
        end%
end%_______________________________________________|
%Ksis:Global system stiffness matrix
clear n sat sut 


disp('Applying Essential Boundary Conditions')

% Apply Essential Boundary Conditions

bcwt=trace(Ksis)/24;

P =P'- Ksis(:,dofConstraints(:,1))*dofConstraints(:,2);
Ksis(:,dofConstraints(:,1)) = 0;
Ksis(dofConstraints(:,1),:) = 0;
Ksis(dofConstraints(:,1),dofConstraints(:,1)) = bcwt*speye(length(dofConstraints(:,1)));
P(dofConstraints(:,1)) = bcwt*dofConstraints(:,2);



%Ksis(Item,Item)=0;
%___________________________SYSTEM STIFFNESS MATRIX
%for n=1:No;
 %       for sat=1:24;
  %          for sut=1:24;
   %            if (R(n,sat)~=0)
    %              if (R(n,sut)~=0);
     %               Ksis(R(n,sut),R(n,sat))=Ksis(R(n,sut),R(n,sat)) + K(sat,sut,n);
      %            end%    
       %        end%
        %    end%
        %end%
%end%_______________________________________________|
%Ksis:Global system stiffness matrix
%clear n sat sut 

disp('Solving')


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
    %Ku = inv(Ksis);
    D = sparseKsis\P;

clear equation Ku

toc

%_________________________________MODAL DISPLACEMENT
 for  s = 1 : No;
        for m = 1 :24;
             u = R(s, m);
             if u ~=0; Hu(s,m,:) = D(u) ;
             else      Hu(s,m,:)  = 0   ;
             end%
        end%   
end%__________________________________________________|
%Hu:per element node`s displacements    
clear s m u 


%Post processing

%____FUNCTION____:Element shape functions partial derivatives.

for s=1:No

    
                    
    
%Moving global cartesian nodes to element nodes
for u=1:8
pxn(u,:) =Cor(Pos(s,u),1);  %x(u)
pyn(u,:) =Cor(Pos(s,u),2);  %y(u)
pzn(u,:) =Cor(Pos(s,u),3);  %z(u)
end

e=0;
n=0;
J=0;


[Jacobi,Hx,Hy,Hz]=AconnectH8(pxn,pyn,pzn,e,n,J);

%[A]=[Equbilibrium differantial matrix].[shell shape matrix]
%[A]=[d].[N]

for f=1:8;
A(:,3*(f-1)+1,s)= [Hx(f);   0   ;   0  ;  Hy(f) ;   0   ;  Hz(f)]; 
A(:,3*(f-1)+2,s)= [ 0   ; Hy(f) ;   0  ;  Hx(f) ;  Hz(f);   0   ]; 
A(:,3*(f-1)+3,s)= [ 0   ;   0   ; Hz(f);   0    ;  Hy(f);  Hx(f)]; 
end

            

end

clear  f u s

for s=1:No
   
        strain(:,s)=A(:,:,s)*Hu(s,:)';
    
       
end
clear s 
%Plot(strain(1,:));

for s=1:No
X(:,s)=Cor(Pos(s,:)',1);
end

for s=1:No
Y(:,s)=Cor(Pos(s,:)',2);
end

for s=1:No
Z(:,s)=Cor(Pos(s,:)',3);
end


    
clear node
xcen=(mean(X));
ycen=mean(Y);
zcen=mean(Z);




%__________________________________________________STRESS AND STRAIN ANALYSIS    
    for s=1:No;                                                            
        for i=1:8;%Element total node number                                
            switch(i)
                    case {1},
                    %       x=-a/2; y=-b/2; z=-c/2;  %Element node for (1)
                             e=-1  ; n=-1 ; J=-1 ;  %Homogen coordinates
                    case {2},
                    %        x=a/2 ; y=-b/2; z=-c/2;  %Element node for (2)
                             e= 1   ; n=-1 ; J=-1 ;  %Homogen coordinates
                    case {3},                    
                    %        x=a/2;y=b/2; z=-c/2;    %Element node for (3)
                             e= 1  ; n= 1 ; J=-1 ;   %Homogen coordinates
                    case {4},
                    %        x=-a/2;y=b/2;z=-c/2;   %Element node for (4)
                             e=-1  ; n= 1 ; J=-1; %Homogen coordinates
                    
                    case {5}, 
                        %    x=-a/2; y=-b/2;z=-c/2;  %Element node for (5)
                             e=-1  ; n=-1 ; J=1;     %Homogen coordinates
                    
                    case {6}, 
                        %    x=a/2 ; y=-b/2; z=c/2;  %Element node for (6)
                             e=1; n=-1; J=1;         %Homogen coordinates
                
                    case {7}, 
                         %   x=a/2;y=b/2; z=c/2     %Element node for (7)
                             e=1; n=1; J=1;         %Homogen coordinates
                        
                    case {8}, 
                        %    x=-a/2;y=b/2; z=c/2   %Element node for (8)
                             e=-1; n=1 ; J=1;      %Homogen coordinates
                       
              end

%        x=gridspacex/2*e;%parametric axis convert homogen axis
%        y=gridspacey/2*n;

%____FUNCTION____:Element shape functions partial derivatives.
[Jacobi,Hx,Hy,Hz]=AconnectH8(pxn,pyn,pzn,e,n,J);


%Element Connection Matrix [A]
for f=1:8;

A(:,3*(f-1)+1,s)= [Hx(f);   0   ;   0  ;  Hy(f) ;   0   ;  Hz(f)]; 
A(:,3*(f-1)+2,s)= [ 0   ; Hy(f) ;   0  ;  Hx(f) ;  Hz(f);   0   ]; 
A(:,3*(f-1)+3,s)= [ 0   ;   0   ; Hz(f);   0    ;  Hy(f);  Hx(f)]; 

end
        
        Qeglobal(:,i,s)=A(:,:,s)*Hu(s,:)';   %Element global node`s Strains
        Qsglobal(:,i,s)=C*Qeglobal(:,i,s);   %Element global node`s Stresses
        
%_____________________________________________Exstrem stress values
        sigmax=Qsglobal(1,i,s);
        sigmay=Qsglobal(2,i,s);
        sigmaz=Qsglobal(3,i,s);
         tauxy=Qsglobal(4,i,s);
         tauyz=Qsglobal(5,i,s);
         tauxz=Qsglobal(6,i,s);       

[taumax,taumin,sigman1,sigman2,sigman3,sigmavonmises]= ...
threedimensionalstress(sigmax,sigmay,sigmaz,tauxy,tauyz,tauxz,1);

%[taumax,taumin,sigman1,sigman2,sigman3,sigmavonmises]= ...
%threedimensionalstress(sigmax,sigmay,sigmaz,tauxy,tauyz,tauxz,0)

Qss(:,i,s)=[taumax
            taumin
            sigman1
            sigman2
            sigman3
            sigmavonmises];
%_____________________________|


        
        end
    end




% Put X, Y, Z co-ordinates into meshgrid form

%for k=1:(meshz+1)
 %   for j=1:(meshy+1)
  %      for i=1:(meshx+1)
   %         X(i,j,k)=Cor(i+(j-1)*(meshx+1)+(k-1)*(meshx+1)*(meshy+1),1);
    %        Y(i,j,k)=Cor(i+(j-1)*(meshx+1)+(k-1)*(meshx+1)*(meshy+1),2);
     %  end
    %end    
%end

% Put strain data into slice form

%node=1;

%############# GRAPHICAL INTERFACE

px=1:meshx+1;%meshgridx 
py=1:meshy+1;%meshgridy 
pz=1:meshz+1;


for n=1:No
for i=1:8 %Element Nodes
    V1(Pos(n,i),:)=Qeglobal(1,i,n);%epsilonxx
    V2(Pos(n,i),:)=Qeglobal(2,i,n);%epsilonyy
    V3(Pos(n,i),:)=Qeglobal(3,i,n);%epsilonzz
    V4(Pos(n,i),:)=Qeglobal(4,i,n);
    V5(Pos(n,i),:)=Qeglobal(5,i,n);
    V6(Pos(n,i),:)=Qeglobal(6,i,n);
end
end


%In the this area element stress values is depend element mesh functions
nm=1;


 for k=1:meshz+1
    for j=1:meshy+1 ;   
        for i=1:meshx+1 ;
   
        
            Qe1(i,j,k,:)=V1(nm) ;
            Qe2(i,j,k,:)=V2(nm) ;
            Qe3(i,j,k,:)=V3(nm) ;
            Qe4(i,j,k,:)=V4(nm) ;
            Qe5(i,j,k,:)=V5(nm) ;
            Qe6(i,j,k,:)=V6(nm) ;
            nm=nm+1;
        
        end
    end
end

dimensions=[totala,totalb,totalc,meshx,meshy,meshz,No];

Qe=[Qe1,Qe2,Qe3,Qe4,Qe5,Qe6];


dlmwrite('dimensions.dat_t2',dimensions);
dlmwrite('strain.dat_t2',strain);
dlmwrite('X.dat_t2',X);
dlmwrite('Y.dat_t2',Y);
dlmwrite('Z.dat_t2',Z);
dlmwrite('xcen.dat_t2',xcen);
dlmwrite('ycen.dat_t2',ycen);
dlmwrite('zcen.dat_t2',zcen);
dlmwrite('Qe.dat_t2',Qe);
dlmwrite('Hx.dat_t2',Hx);
dlmwrite('Hy.dat_t2',Hy);
dlmwrite('Hz.dat_t2',Hz);
dlmwrite('C.dat_t2',C);
dlmwrite('Cor.dat_t2',Cor);
dlmwrite('Re.dat_t2',Re);
dlmwrite('R.dat_t2',R);
dlmwrite('Qeglobal.dat_t2',Qeglobal);
dlmwrite('Qsglobal.dat_t2',Qsglobal);


disp('Programme Finished')








