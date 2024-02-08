function [x, y] = DMRGIdx2Coor(Ly, Lx, idx)
x_grain = floor(idx/(3*Ly-1));
x = 2 *x_grain;
if( mod(idx, 3*Ly-1) >= 2*Ly-1)
    x = x+1;
    y = 2*(mod(idx, 3*Ly-1) - (2*Ly-1));
else
    y =  mod(idx, 3 *Ly-1);
end
end