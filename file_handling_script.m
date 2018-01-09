% Construct Path for data file corresponding to these parameters

filepath=Filepath(handles.nelx,handles.nely,handles.nelz,handles.aspectab,handles.aspectac,handles.bcmode,handles.ex,handles.ey,handles.ez,handles.styecast_ratio);
filename=Filename(handles.nelx,handles.nely,handles.nelz,handles.aspectab,handles.aspectac,handles.bcmode,handles.ex,handles.ey,handles.ez,handles.styecast_ratio);


% Now check whether this data file allready exists

file_existence=Does_file_exist(filepath);

% Based on this information we can now decide either read in the file or
% call the FEA routine to get the data needed to create the file

if file_existence==true
    
    load(eval(['' 'filepath' '']));
    
    SX=eval([filename '.SX']);
    SY=eval([filename '.SY']);
    SZ=eval([filename '.SZ']);
    SSTRAIN=eval([filename '.SSTRAIN']);
    X=eval([filename '.X']);
    Y=eval([filename '.Y']);
    Z=eval([filename '.Z']);
    Disp=eval([filename '.Disp']);
    dimensions=eval([filename '.dimensions']);
    
    
else
    
    [ SX,SY,SZ,SSTRAIN,X,Y,Z,Disp,dimensions]=FEA_Routine(handles.nelx,handles.nely,handles.nelz,handles.aspectab,handles.aspectac,handles.bcmode,handles.ex,handles.ey,handles.ez,2,handles.styecast_ratio);
 
    eval([filename '=struct( ''' 'SX' ''' ,SX, ''' 'SY' ''' ,SY,''' 'SZ' ''',SZ, ''' 'SSTRAIN' ''' ,SSTRAIN, ''' 'X' ''' ,X, ''' 'Y' ''' ,Y, ''' 'Z' ''' ,Z, ''' 'Disp' ''' ,Disp, ''' 'dimensions' ''' ,dimensions);']);
    
   
 
    save(filepath,filename);
 
   
end