function [x_idx, y_idx, x_coor, y_coor] = PEPSIdx2Coor(Lx, idx)
% triangle lattice
y_idx = fix(idx/Lx);
x_idx = mod(idx, Lx);

a = 1; %lattice constant
x_coor = x_idx * a + y_idx * a/2;
y_coor = y_idx * a * sqrt(3)/2;
end