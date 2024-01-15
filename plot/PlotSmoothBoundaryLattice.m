function [] = PlotSmoothBoundaryLattice(Ly, Lx, plot_bond, plot_circle)
if(nargin < 3)
    plot_bond = 1;
end
if(nargin < 4)
    plot_circle =0;
end
N = 3 * Lx * Ly - Lx - Ly;


a = 1; % lattice constant

lattice_linewidth = 3;
site_size = 80;
bond_color = [233, 196, 107]/256;
site_color = [042, 157, 142]/256;
circle_color = [164, 005, 069]/256;
each_row_site_num = 3*Lx-1;

[x_set, y_set] = GenerateSmoothBoundaryLatticeCoord(Ly, Lx);


if(plot_bond)
    for site_idx = 0:N-1
        y = floor(site_idx / each_row_site_num);
        x_in_row = mod(site_idx, each_row_site_num);
        if(y < Ly -1)
            x_grain = floor(x_in_row /3);
            x_inner = mod(x_in_row, 3);
            if(x_inner == 0)
                if(x_grain < Lx -1)
                    line([x_set(site_idx+1),x_set(site_idx+1)+a],[y_set(site_idx+1),y_set(site_idx+1)],'color',bond_color,'linewidth',lattice_linewidth);  hold on;
                end
                line([x_set(site_idx+1),x_set(site_idx+1)+a/2],[y_set(site_idx+1),y_set(site_idx+1)+sqrt(3)/2 *a ],'color',bond_color,'linewidth',lattice_linewidth);  hold on;
            elseif(x_inner==1)
                line([x_set(site_idx+1),x_set(site_idx+1)+a/2],[y_set(site_idx+1),y_set(site_idx+1)+sqrt(3)/2 *a ],'color',bond_color,'linewidth',lattice_linewidth);  hold on;
                if(x_grain > 0)
                    line([x_set(site_idx+1),x_set(site_idx+1)-a/2],[y_set(site_idx+1),y_set(site_idx+1)+sqrt(3)/2 *a ],'color',bond_color,'linewidth',lattice_linewidth);  hold on;
                end
            else
                line([x_set(site_idx+1),x_set(site_idx+1)+a],[y_set(site_idx+1),y_set(site_idx+1)],'color',bond_color,'linewidth',lattice_linewidth);  hold on;
                line([x_set(site_idx+1),x_set(site_idx+1)-a/2],[y_set(site_idx+1),y_set(site_idx+1)+sqrt(3)/2 *a],'color',bond_color,'linewidth',lattice_linewidth);  hold on;
            end
        else
            if(x_in_row /2  < Lx-1)
                line([x_set(site_idx+1),x_set(site_idx+1)+a],[y_set(site_idx+1),y_set(site_idx+1)],'color',bond_color,'linewidth',lattice_linewidth);  hold on;
            end
        end
    end
end
scatter(x_set,y_set,site_size,site_color,'o','filled');

if(plot_circle)
    for x = 0:Lx-1
        for y = 0:Ly-1
            base_point = [2*a * x+ y, y * sqrt(3)];
            circel_cent = base_point+[a/2, a/sqrt(3)/2];
            viscircles(circel_cent, a/1.2,'Color',circle_color,'LineWidth', 3);
        end
    end
end

axis off;
% axis tight;
axis equal;

set(gca, 'YDir','reverse')

end