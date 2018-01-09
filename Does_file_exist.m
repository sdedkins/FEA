function [ file_existence] = Does_file_exist( filepath)
%Does_file_exist function to determine whether a data file correpsonding to
%a given set of simulation parameters already exists



if exist(filepath, 'file') 
    
    file_existence=true;
else
    file_existence=false;
    
end


end

