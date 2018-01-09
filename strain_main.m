%_________________________________________________________________
%| 8-Node Hexahedral Finite Element Solver
%|
%|
%|
%|
%|
%|
%|
%|
%|
%|
%|
%|
%|                                                              
%|                                                              
%________________________________________________________________
clear
tic
aspectab=10;
aspectac=10;

meshx=100;
meshy=14;
meshz=14;

bcmode=3;

styecastratio=2;

ex=0;
ey=6.62e-4;
ez=1.83e-3;

material_type=1;

nel=meshx*meshy*meshz;    % Number of elements 
nodex=meshx+1;            % Number of Nodes in x-direction
nodey=meshy+1;            % Number of Nodes in y-direction
nodez=meshz+1;            % Number of Nodes in z-direction
nodes=nodex*nodey*nodez;   % Number of Nodes
ndof=3*nodes;             % Number of degrees of freedom


if material_type==1;
    
        E11=2.32;
        E33=2.08;
        E23=0.71;
        E12=1.06;
        E44=0.657;
        E66=0.612;

        E=zeros(6,6);

        E(1,:)=[E11,E12,E23,0,0,0];
        E(2,:)=[E12,E11,E23,0,0,0];
        E(3,:)=[E23,E23,E33,0,0,0];
        E(4,:)=[0,0,0,E44,0,0];
        E(5,:)=[0,0,0,0,E44,0];
        E(6,:)=[0,0,0,0,0,E66];

        E=E.*1e11;
        
elseif material_type==2
    
        v=0.3;
        Emag=1e11;
    
        E=Emag/((1+v)*(1-2*v))*[1-v ,  v , v ,  0      ,    0    ,    0     ;   
                              v  ,1-v , v ,  0      ,    0    ,    0     ;  
                              v  ,  v ,1-v,  0      ,    0    ,    0     ; 
                              0  ,  0 , 0 ,(1-2*v)/2,    0    ,    0     ;
                              0  ,  0 , 0 ,  0      ,(1-2*v)/2,    0     ;
                              0  ,  0 , 0 ,  0      ,    0    ,(1-2*v)/2];
                  
else
    
    error('Unrecognised value for material type');

end


[Cor , Pos, sample,total]=Preprocessor( aspectab,aspectac,meshx,meshy,meshz,bcmode,styecastratio,nodes,nodex,nodey,nodez);

[Re,R,dofConstraints,P]=Boundary_Conditions(Pos,Cor,ex, ey, ez,nodes,nel,ndof,bcmode,sample,total);

[D,sparseKsis]=Solver_v2(Cor,Pos,R,E,dofConstraints,P,nel,ndof);

 [strain]=Strain_calc(Cor,Pos,R,nel,D);        
 
 
 
 for s=1:nel
X(:,s)=Cor(Pos(s,:)',1);
end

for s=1:nel
Y(:,s)=Cor(Pos(s,:)',2);
end

for s=1:nel
Z(:,s)=Cor(Pos(s,:)',3);
end

xcen=(mean(X));
ycen=mean(Y);
zcen=mean(Z);

dimensions=[total(1),total(2),total(3),meshx,meshy,meshz,nel];

dlmwrite('dimensions.dat_t2',dimensions);
dlmwrite('strain.dat_t2',strain);
dlmwrite('X.dat_t2',X);
dlmwrite('Y.dat_t2',Y);
dlmwrite('Z.dat_t2',Z);
dlmwrite('xcen.dat_t2',xcen);
dlmwrite('ycen.dat_t2',ycen);
dlmwrite('zcen.dat_t2',zcen);
dlmwrite('Cor.dat_t2',Cor);
dlmwrite('Re.dat_t2',Re);
dlmwrite('R.dat_t2',R);



toc
disp('Programme Finished')
        
        
        

               
               



