Ly = 4;
Lx = 12;
mark_energy_value = 0;
energy_text_size=10;

bond_with_scale = 10;
postive_bond_color = [233, 196, 107]/256;
minus_bond_color = [042, 157, 142]/256;

bond_energy_file1 = ['../data/kagome_dmrg', num2str(Ly), 'x', num2str(Lx),'bondzz.json'];
bond_energy_data1 = jsondecode(fileread(bond_energy_file1));
bond_energy_file2 = ['../data/kagome_dmrg', num2str(Ly), 'x', num2str(Lx),'bondpm.json'];
bond_energy_data2 = jsondecode(fileread(bond_energy_file2));
bond_energy_file3 = ['../data/kagome_dmrg', num2str(Ly), 'x', num2str(Lx),'bondmp.json'];
bond_energy_data3 = jsondecode(fileread(bond_energy_file3));

bond_num = numel(bond_energy_data1);

% PlotSmoothBoundaryLattice(Ly, Lx, 1);
a=1;
for i = 1:bond_num
    site1_idx = bond_energy_data1{i}{1}(1);
    site2_idx = bond_energy_data1{i}{1}(2);
    bond_energy = bond_energy_data1{i}{2} + 1/2 * bond_energy_data2{i}{2} + 1/2 * bond_energy_data3{i}{2};

    [x1_coor, y1_coor] = DMRGIdx2Coor(Ly, Lx, site1_idx);
    [x2_coor, y2_coor] = DMRGIdx2Coor(Ly, Lx, site2_idx);

    x1_pos = x1_coor * a + y1_coor/2 * a;
    y1_pos = y1_coor * sqrt(3)/2 * a;
    x2_pos = x2_coor * a + y2_coor/2 * a;
    y2_pos = y2_coor * sqrt(3)/2 * a;
    if(bond_energy > 0)
        line([x1_pos,x2_pos],[y1_pos,y2_pos],'color',postive_bond_color,'linewidth',bond_with_scale* bond_energy);  hold on;
    else
        line([x1_pos,x2_pos],[y1_pos,y2_pos],'color',minus_bond_color,'linewidth',-bond_energy* bond_with_scale);  hold on;
    end

    if(mark_energy_value)
        T = text((x1_pos+x2_pos)/2, (y1_pos + y2_pos)/2, [num2str(bond_energy)],'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

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
