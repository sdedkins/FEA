function [ sx,sy,sz ] = Outside_Surface( SX,SY,SZ )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%% Define outside planes of sample for plotting____%%   
sx=[SX(1,1,1), SX(size(SX,1),size(SX,2),size(SX,3))];
sy=[SY(1,1,1),SY(size(SY,1),1,1)];
sz=[SZ(1,1,1), SZ(size(SZ,1),size(SZ,2),size(SZ,3))];



end

