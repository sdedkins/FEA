clear 

close all hidden


dimensions=dlmread('dimensions.dat_t2');
strain=dlmread('strain.dat_t2');
X=dlmread('X.dat_t2');
Y=dlmread('Y.dat_t2');
Z=dlmread('Z.dat_t2');
xcen=dlmread('xcen.dat_t2');
ycen=dlmread('ycen.dat_t2');
zcen=dlmread('zcen.dat_t2');


[SX,SY,SZ,sx,sy,sz,SSTRAIN]=Sliceform(dimensions,strain,xcen,ycen,zcen);

[clinex, clines]=Cline(SX,SY,SZ,SSTRAIN);

%% Plot using slice function____%%

figure(1) 


slice(SX,SY,SZ,SSTRAIN(:,:,:,2),sx,sy,sz);
shading interp;
axis equal;
axis vis3d;
colormap jet;
colorbar;
xlabel('x/mm');
ylabel('y/mm');
zlabel('z/mm');

%%_____Plot Strains along the centreline_____%%

figure(2) 

subplot(1,3,1), plot(clinex./100,clines(:,1)), axis square, xlabel('x/mm'), ylabel('\epsilon_{xx}');
subplot(1,3,2), plot(clinex./100,clines(:,2)), axis square ,xlabel('x/mm'), ylabel('\epsilon_{yy}');
subplot(1,3,3), plot(clinex./100,clines(:,3)), axis square, xlabel('x/mm'), ylabel('\epsilon_{zz}');


%
figure(3) 

subplot

subplot(1,3,1), plot(clinex./100,clines(:,4)./2), axis square, xlabel('x/mm'), ylabel('\epsilon_{xy}');
subplot(1,3,2), plot(clinex./100,clines(:,5)./2), axis square ,xlabel('x/mm'), ylabel('\epsilon_{yz}');
subplot(1,3,3), plot(clinex./100,clines(:,6)./2), axis square, xlabel('x/mm'), ylabel('\epsilon_{zx}');



