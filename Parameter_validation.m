function [ validated ] = Parameter_validation( meshx,meshy,meshz,aspectab,aspectac,bcmode,styecast_ratio)
%Parameter_validation: Returns 0 if parameters are valid and 0 if not.
%   Detailed explanation goes here

validated=0;

if mod(meshx,1)~=0
    validated=1;
    display('Non-Integer Element Number Entered')
end


if mod(meshy,1)~=0
    validated=1;
    display('Non-Integer Element Number Entered')
end


if mod(meshz,1)~=0
    validated=1;
    display('Non-Integer Element Number Entered')
end

if aspectab<=0
    validated=1;
    display('Negative Aspect Ratio Entered')
end

if aspectac<=0
    validated=1;
    display('Negative Aspect Ratio Entered')
end

if mod(styecast_ratio,1) ~=0
    validated=1;
    display('Non-Integer Styecast Ratio Entered')
end

if (bcmode==1 || bcmode==2 || bcmode==3)
    
else
    
    display('Unrecognised Boundary Conditions')
    validated=1;
end


    


end

