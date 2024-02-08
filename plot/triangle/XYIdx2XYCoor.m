function [x_coor, y_coor] = XYIdx2XYCoor(x_idx, y_idx, a)
%XYIDX2XYCOOR Summary of this function goes here
%   Detailed explanation goes here
x_coor = x_idx * a + y_idx * a /2;
y_coor = y_idx * a * sqrt(3)/2;
end

