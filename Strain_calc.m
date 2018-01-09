function [ strain,Hu ] = Strain_calc(Cor,Pos,R,nel,D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



%_________________________________MODAL DISPLACEMENT
 for  s = 1 : nel;
        for m = 1 :24;
             
             Hu(s,m,:) = D(R(s,m)) ;
            
        end
        
end%__________________________________________________|
%Hu:per element node`s displacements    
clear s m u 



%____FUNCTION____:Element shape functions partial derivatives.

for s=1:nel

    
                    
    
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

for s=1:nel
   
        strain(:,s)=A(:,:,s)*Hu(s,:)';
    
       
end
clear s 




end

