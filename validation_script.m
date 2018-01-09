validated=Parameter_validation( handles.nelx,handles.nely,handles.nelz,handles.aspectab,handles.aspectac,handles.bcmode,handles.styecast_ratio);

% Tell user paramters validated or give error message

if validated==0
    
    display ('Parameters Validated')
    
else
    
    error('Parameters could not be validated');
    
end

clear validated;


