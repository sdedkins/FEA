clear

close all hidden

% Vector holding which aspect ratios data should be plotted for
aspect=[10,9,8,7,6,5,4,3,2.7,2.5];

%Simulation Parameters
meshx_sample=60;
meshy=15;
meshz=15;
thermal_mismatch=-0.2;
applied_strain=0.05;

% Call function to run simulation or read in data file for each set of 
%simulation conditions and return a 3-D structure array of the form 
%H_Data(i, bcmode, matched)= structure holding strain data for aspect
%ratio i, boundary condition mode bcmode, and matched=1 -> thermal
%mismatch matched=2 -> thermally matched

[filenames,filepaths,H_DATA]=Homog_Get(aspect,meshx_sample,meshy,meshz,...
    thermal_mismatch,applied_strain);

% Call to function that extracts the co-ordinates diagonal strains at the...
%centre of the sample in each simulation and returns as a structure array
%of the form CENTRE_STRAIN(i,bcmode,matched)= structure corresponding to
%those simulation parameters. Each structure has the form
%struct('x',{},'y',{},'z',{},'exx',{},'eyy',{},'ezz',{})

CENTRE_STRAIN=Centre_strain(H_DATA,aspect);




for i=1:length(aspect)
    for bcmode=1:3
        for matched=1:2
            
    exx(i,bcmode,matched)=CENTRE_STRAIN(i,bcmode,matched).exx*100;
    eyy(i,bcmode,matched)=CENTRE_STRAIN(i,bcmode,matched).eyy*100;
    ezz(i,bcmode,matched)=CENTRE_STRAIN(i,bcmode,matched).ezz*100;
    
        end
    end
end


figure
subplot(2,2,1);
plot(aspect,exx(:,2,1),'o',aspect,eyy(:,2,1),'+',aspect,ezz(:,2,1),'x');
xlabel('Aspect Ratio');
ylabel('\epsilon /%');
title({'Strain at Centre','Ends is stycast with 0.2% thermal mismatch.'});
legend('\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}','location','southeast');

subplot(2,2,2);
plot(aspect,-eyy(:,2,1)./exx(:,2,1),'o');
ylabel('$$\nu_{eff}$$','interpreter','latex');
title({'Effective Poisson Ratio.','Ends is stycast with 0.2% thermal mismatch.'}); 

subplot(2,2,3);
plot(aspect,exx(:,2,2)/applied_strain,'o',aspect,eyy(:,2,2)/(-0.3*applied_strain),'+',aspect,ezz(:,2,2)/-(0.3*applied_strain),'x');
xlabel('Aspect Ratio');
ylabel('$$\frac{\epsilon}{\epsilon_{ideal}}$$','interpreter','latex');
title({'Ratio of actual strain at centre to strain for ideal uniaxial case.','Ends is stycast with no thermal mismatch.'});
legend('\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}','location','southeast');


subplot(2,2,4);
plot(aspect,(applied_strain./exx(:,2,2)-1)*100,'o',aspect,((-0.3*applied_strain)./eyy(:,2,2)-1)*100,'+',aspect,(-(0.3*applied_strain)./ezz(:,2,2)-1)*100,'x');
xlabel('Aspect Ratio');
ylabel('$$\frac{L_{eff}-L}{L} / \%$$','interpreter','latex');
title({'Effective length of sample.','Ends is stycast with no thermal mismatch.'}); 
legend('\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}','location','northeast');


