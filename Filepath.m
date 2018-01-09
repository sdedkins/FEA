function [ filepath ] = Filepath( meshx,meshy,meshz,aspectab,aspectac,bcmode,ex,ey,ez,styecast_ratio )
%Filepath: Returns file path for data file where data for a given set of
%simulation parameters is stored
%   
%Returns file path in the form of a string
%'\Data\meshx_meshy_meshz_aspectab,_aspectac_ab

str_meshx=num2str(meshx, '%u');
str_meshy=num2str(meshy, '%u');
str_meshz=num2str(meshz, '%u');

str_aspectab=num2str(aspectab, '%tx');
str_aspectac=num2str(aspectac, '%tx');

str_bcmode=num2str(bcmode, '%u');

str_styecast_ratio=num2str(styecast_ratio, '%u');

str_ex=num2str(ex, '%tx');
str_ey=num2str(ey, '%tx');
str_ez=num2str(ez, '%tx');


filepath=strcat('Data\FEA',str_meshx,'_',str_meshy,'_',str_meshz,'_',str_aspectab,'_',str_aspectac,'_',str_styecast_ratio,'_',str_bcmode,'_',str_ex,'_',str_ey,'_',str_ez,'.mat');

end

