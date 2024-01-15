function [] = PlotBondEnergy( energy_bond_info, bond_energys, mark_energy_value)
% energy_bond_info :[x1, y1, x2,y2]
if(nargin < 3)
    mark_energy_value = 0;
end
bond_num = size(energy_bond_info, 1);
Ly = (max(energy_bond_info(:,1))+2)/2;
Lx = (max(energy_bond_info(:,2))+2)/2;

% PlotSmoothBoundaryLattice(Ly, Lx, 0);
energy_text_size = 40;

a=1;

bond_with_scale = 10;
postive_bond_color = [233, 196, 107]/256;
minus_bond_color = [042, 157, 142]/256;

for i = 1:bond_num
    y1_coor = energy_bond_info(i, 1);
    x1_coor = energy_bond_info(i, 2);
    y2_coor = energy_bond_info(i, 3);
    x2_coor = energy_bond_info(i, 4);

    x1_pos = x1_coor * a + y1_coor/2 * a;
    y1_pos = y1_coor * sqrt(3)/2 * a;
    x2_pos = x2_coor * a + y2_coor/2 * a;
    y2_pos = y2_coor * sqrt(3)/2 * a;
    if(bond_energys(i) > 0)
        line([x1_pos,x2_pos],[y1_pos,y2_pos],'color',postive_bond_color,'linewidth',bond_with_scale* bond_energys(i));  hold on;
    else
        line([x1_pos,x2_pos],[y1_pos,y2_pos],'color',minus_bond_color,'linewidth',-bond_energys(i)* bond_with_scale);  hold on;
    end
    if(mark_energy_value)
        T = text((x1_pos+x2_pos)/2, (y1_pos + y2_pos)/2, [num2str(bond_energys(i))],'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

        set(T,'Interpreter','latex');
        set(T,'Fontsize',energy_text_size);

        if(x1_coor == x2_coor && y2_coor == y1_coor+1)
            set(T,'Rotation',-60);
        elseif(x1_coor+1 == x2_coor && y2_coor == y1_coor-1)
            set(T,'Rotation',60);
        end
    end
end

axis off;
% axis tight;
axis equal;

set(gca, 'YDir','reverse')

end