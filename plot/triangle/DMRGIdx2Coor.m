function [x_idx, y_idx, x_coor, y_coor] = DMRGIdx2Coor(Ly, idx)
% triangle lattice
x_idx = fix(idx/Ly);
y_idx = mod(idx, Ly);

a = 1; %lattice constant
x_coor = x_idx * a + y_idx * a/2;
y_coor = y_idx * a * sqrt(3)/2;
end