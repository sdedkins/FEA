function [ filename,filepath,H_DATA ] = Homog_Get( aspect,meshx_sample,...
    meshy,meshz,thermal_mismatch,applied_strain)
%H Summary of this function goes here
%   Detailed explanation goes here


no_aspect=length(aspect);
stycast_ratio=2;

for i=1:no_aspect
        for bcmode=1:3
            for matched=1:2  % 1 corresponds to thermal mismatch, 
                %2 corresponds to thermally matched 

%_____Calculate Required Strains________________________________________%

ex=applied_strain/100;

if matched==1 
    
    ex=0;
   
    if bcmode==1||bcmode==2
    
        
            ey=thermal_mismatch/100;
            ez=thermal_mismatch/100;
            
    elseif bcmode==3
        
        ey=thermal_mismatch/100;
        ez=0;
        
    else error('Unrecognised Boundary Condition in Strain Calculation')
    end
    
elseif matched==2
    
    ex=applied_strain/100;
    
    if bcmode==1||bcmode==2
        
        ey=0;
        ez=0;
    elseif bcmode==3
        
        ey=0;
        ez=0;
        
    else error('Unrecognised Boundary Condition in Strain Calculation')
    end
        
else error('The variable <<matched>> has an unexpected value')  
      
            
            
    
    end
    

                
                

%_____________________ Parameter Validation____________________________

validated=Parameter_validation( meshx_sample,meshy,meshz,aspect(i),...
    aspect(i),bcmode,stycast_ratio);

% Tell user paramters validated or give error message

if validated==0
    
    display ('Parameters Validated')
    
else
    
    error('Parameters could not be validated');
    
end

clear validated;

%________________________________________________________________________

%____________________File Handling_______________________________________

% Construct Path for data file corresponding to these parameters

filepath(i,bcmode,matched)=cellstr(Filepath(meshx_sample,meshy,meshz,...
    aspect(i),aspect(i),bcmode,ex,ey,ez,stycast_ratio));
filename(i,bcmode,matched)=cellstr(Filename(meshx_sample,meshy,...
    meshz,aspect(i),aspect(i),bcmode,ex,ey,ez,stycast_ratio));


% Now check whether this data file allready exists

file_existence=Does_file_exist(filepath{i,bcmode,matched});

% Based on this information we can now decide to either read in the file or
% call the FEA routine to get the data needed to create the file

if file_existence==true
    
    load(eval(['' 'filepath{i,bcmode,matched}' '']));
    
    SX=eval([filename{i,bcmode,matched} '.SX']);
    SY=eval([filename{i,bcmode,matched} '.SY']);
    SZ=eval([filename{i,bcmode,matched} '.SZ']);
    SSTRAIN=eval([filename{i,bcmode,matched} '.SSTRAIN']);
    X=eval([filename{i,bcmode,matched} '.X']);
    Y=eval([filename{i,bcmode,matched} '.Y']);
    Z=eval([filename{i,bcmode,matched} '.Z']);
    Disp=eval([filename{i,bcmode,matched} '.Disp']);
    dimensions=eval([filename{i,bcmode,matched} '.dimensions']);
    
    % Write structure containing data to a structure array for return
    H_DATA(i,bcmode,matched)=eval(filename{i,bcmode,matched});
    
else
    
    %Call FEA Routine
    [ SX,SY,SZ,SSTRAIN,X,Y,Z,Disp,dimensions]=FEA_Routine(meshx_sample,...
        meshy,meshz,aspect(i),aspect(i),bcmode,ex,ey,ez,2,stycast_ratio);
 %Save data in structure with name of intended file
    eval([filename{i,bcmode,matched} '=struct( ''' 'SX' ''' ,SX, ''' 'SY' ''' ,SY,''' 'SZ' ''',SZ, ''' 'SSTRAIN' ''' ,SSTRAIN, ''' 'X' ''' ,X, ''' 'Y' ''' ,Y, ''' 'Z' ''' ,Z, ''' 'Disp' ''' ,Disp, ''' 'dimensions' ''' ,dimensions);']);
    
    %Save structure to a .mat file
    save(filepath{i,bcmode,matched},filename{i,bcmode,matched});
    % Write structure containing data to a structure array for return
    H_DATA(i,bcmode,matched)=eval(filename{i,bcmode,matched});

end

            end
        end     
end
       
end

%________________________________________________________________________%




