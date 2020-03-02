function y = from1d2d(input,fullsize)
%from1d2d Calculate 2D rotationally symmetrical image extended from 1D profile
% SYNOPSIS:
% y = from1d2d(input,fullsize)
%
% PARAMETERS:
%    input: required input is vector or image with size 1xn
% fullsize: 
%
% OUTPUT:
%   y: 2D image

ra = rr(fullsize,fullsize);
maxr = ceil(max(ra));
cut = fullsize/2-length(input);
t1 = 1:length(input);
out1d_ext = [zeros(1,cut),input,zeros(1,maxr+1-max(t1))];
% Calculate 2D PSD extended from 1D average
y=lut(ra,out1d_ext);
