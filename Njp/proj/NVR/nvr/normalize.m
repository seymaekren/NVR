function [n] = normalize(v)
% This function accepts a vector v and normalizes it to length one
%
n = v / (sqrt(sum (v.*v)));
