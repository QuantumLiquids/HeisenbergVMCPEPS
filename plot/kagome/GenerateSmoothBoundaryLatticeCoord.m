function [x_set, y_set] = GenerateSmoothBoundaryLatticeCoord(Ly, Lx)
N = 3 * Lx * Ly - Lx - Ly;


x_set = zeros(1,N);
y_set = zeros(1,N);

a = 1; % lattice constant

each_row_site_num = 3*Lx-1;
for site_idx = 0:N-1
    % row major
    y = floor(site_idx / each_row_site_num);
    x_in_row = mod(site_idx, each_row_site_num);
    if(y < Ly -1)
        x_grain = floor(x_in_row /3);
        x_inner = mod(x_in_row, 3);
        if(x_inner == 0)
            x_set(site_idx+1) = x_grain * (2*a) + y * a;
            y_set(site_idx+1) = y * (sqrt(3) * a);
        elseif(x_inner==1)
            x_set(site_idx+1) = x_grain * (2*a) + y * a + a/2;
            y_set(site_idx+1) = y * (sqrt(3) * a) + sqrt(3)/2 * a;
        else
            x_set(site_idx+1) = x_grain * (2*a) + y * a+a;
            y_set(site_idx+1) = y * (sqrt(3) * a);
        end
    else
        x_set(site_idx+1) = x_in_row * a + y * a;
        y_set(site_idx+1) = y * (sqrt(3) * a);
    end
end


end