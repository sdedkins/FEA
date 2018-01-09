clear
clf

dimensions=dlmread('dimensions.dat_t2');
strain=dlmread('strain.dat_t2');
X=dlmread('X.dat_t2');
Y=dlmread('Y.dat_t2');
Z=dlmread('Z.dat_t2');
xcen=dlmread('xcen.dat_t2');
ycen=dlmread('ycen.dat_t2');
zcen=dlmread('zcen.dat_t2');

%[totala,totalb,totalc,meshx,meshy,meshz,No]= dimensions(1,:);
meshx=dimensions(4);
meshy=dimensions(5);
meshz=dimensions(6);


clinex=(xcen(2201:2250)+xcen(2251:2300))./2;
cliney=(ycen(2201:2250)+ycen(2251:2300))./2;
clines=(strain(:,2201:2250)+strain(:,2251:2300)+strain(:,2701:2750)+strain(:,2751:2800))./4;

% plot strain along the centre line of the block

L=polarmap;

%figure 

%colormap(L);

%subplot(3,1,1), fill(X(1:4,1:(meshx*meshy))./100,Y(1:4,1:(meshx*meshy))./100,strain(4,1:(meshx*meshy))), axis equal, colorbar, xlabel('x/mm'), ylabel('y/mm'), zlabel('\epsilon_{xx}');
%subplot(3,1,2), fill(X(1:4,2001:2500)./100,Y(1:4,2001:2500)./100,strain(4,2001:2500)), axis equal , colorbar, xlabel('x/mm'), ylabel('y/mm'), zlabel('\epsilon_{xx}');
%subplot(3,1,3), fill(X(1:4,4501:end)./100,Y(1:4,4501:end)./100,strain(4,4501:end)), axis equal, colorbar, xlabel('x/mm'), ylabel('y/mm'), zlabel('\epsilon_{xx}');


%figure 

%subplot

%subplot(1,3,1), plot(clinex./100,clines(4,:)./2), axis square, xlabel('x/mm'), ylabel('\epsilon_{xy}');
%subplot(1,3,2), plot(clinex./100,clines(5,:)./2), axis square ,xlabel('x/mm'), ylabel('\epsilon_{yz}');
%subplot(1,3,3), plot(clinex./100,clines(6,:)./2), axis square, xlabel('x/mm'), ylabel('\epsilon_{zx}');



figure 

subplot(1,3,1), plot(clinex./100,clines(1,:)), axis square, xlabel('x/mm'), ylabel('\epsilon_{xx}');
subplot(1,3,2), plot(clinex./100,clines(2,:)), axis square ,xlabel('x/mm'), ylabel('\epsilon_{yy}');
subplot(1,3,3), plot(clinex./100,clines(3,:)), axis square, xlabel('x/mm'), ylabel('\epsilon_{zz}');

%figure

%colormap(autumn);

%subplot(3,1,1), fill(Y([1 5 8 4],1:50:4951)./100,Z([1 5 8 4],1:50:4951)./100,strain(1,1:50:4951)), axis equal, colorbar, xlabel('y/mm'), ylabel('z/mm'), zlabel('\epsilon_{xx}');
%subplot(3,1,2), fill(Y([1 5 8 4],1:50:4951)./100,Z([1 5 8 4],1:50:4951)./100,strain(2,1:50:4951)), axis equal, colorbar, xlabel('y/mm'), ylabel('z/mm'), zlabel('\epsilon_{yy}');
%subplot(3,1,3), fill(Y([1 5 8 4],1:50:4951)./100,Z([1 5 8 4],1:50:4951)./100,strain(3,1:50:4951)), axis equal, colorbar, xlabel('y/mm'), ylabel('z/mm'), zlabel('\epsilon_{zz}');


%figure(fig)
%subplot(1,3,1);pcolor(py,px,Qe1(:,:,3,1));title('exx');shading interp; colorbar
%subplot(1,3,2);pcolor(py,px,Qe2(:,:,3,1));title('eyy');shading interp; colorbar;
%subplot(1,3,3);pcolor(py,px,Qe3(:,:,3,1));title('ezz');shading interp; colorbar;



