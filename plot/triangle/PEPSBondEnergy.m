Ly = 10;
Lx = 10;
J2 = 0.125;
Dpeps = 8;
Db = 32;

mark_energy_value  = 0;
energy_text_size = 40;
bond_with_scale = 10;
postive_bond_color = [233, 196, 107]/256;
minus_bond_color = [042, 157, 142]/256;

nn_bond_num =  Lx * (Ly-1)+ (Lx-1) * Ly + (Lx-1) * (Ly-1);
if(J2==0)
    bond_num = nn_bond_num;
else
    bond_num =  Lx * (Ly-1)+ (Lx-1) * Ly + (Lx-1) * (Ly-1) *2 + (Lx-2) * (Ly-1) + (Ly-2)*(Lx-1);
end
site_num = Ly * Lx ;
data_file_name = ['../../data/triangle_energy_statistics', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)];
file_id = fopen(data_file_name,'rb');
energy = fread(file_id, 1, 'double');
en_std = fread(file_id, 1, 'double');
bond_energys = fread(file_id, bond_num, 'double');
fclose(file_id);
fprintf('Energy +- en_std: %f pm %f\n', energy, en_std);

% === Generate the NN bond info ==== %
nn_bond_coor_info = zeros(nn_bond_num, 4); % x1_coor, y1_coor, x2_coor, y2_coor
a = 1; %lattice constant
count = 1;
for y_idx = 0:Ly -1
    % horizontal bond
    for x_idx = 0:Lx-2 %Cpp convention
        [x1_coor, y1_coor] =XYIdx2XYCoor(x_idx, y_idx, a);
        [x2_coor, y2_coor] =XYIdx2XYCoor(x_idx+1, y_idx, a);
        nn_bond_coor_info(count, :)=[x1_coor, y1_coor, x2_coor, y2_coor];
        count = count + 1;
    end
    % diagonal bond
    if(y_idx < Ly - 1)
        for x_idx = 0:Lx-2
            [x1_coor, y1_coor] =XYIdx2XYCoor(x_idx, y_idx+1, a);
            [x2_coor, y2_coor] =XYIdx2XYCoor(x_idx+1, y_idx, a);
            nn_bond_coor_info(count, :)=[x1_coor, y1_coor, x2_coor, y2_coor];
            count = count + 1;
        end
    end
end
for x_idx = 0:Lx-1
    for y_idx = 0:Ly -2
        % vertical bond
        [x1_coor, y1_coor] =XYIdx2XYCoor(x_idx, y_idx, a);
        [x2_coor, y2_coor] =XYIdx2XYCoor(x_idx, y_idx+1, a);
        nn_bond_coor_info(count, :)=[x1_coor, y1_coor, x2_coor, y2_coor];
        count = count + 1;
    end
end

% ==== NN bond Idx Map =====%
% map from all bond energy data to NN bond energy data
nn_bond_idx_map = 1:nn_bond_num;
count1 = 1; % idx to nn bond data
count2 = 1; % idx to row data
if(J2 ~=0)
    for y_idx = 0:Ly -1
        % horizontal bond
        for x_idx = 0:Lx-2 %Cpp convention
            nn_bond_idx_map(count1) = count2;
            count1 = count1 + 1;
            count2 = count2 + 1;
        end

        % diagonal bond
        if(y_idx < Ly - 1)
            for x_idx = 0:Lx-2
                nn_bond_idx_map(count1) = count2;
                count1 = count1 + 1;
                if x_idx < Lx-2
                    count2 = count2 + 3;
                else
                    count2 = count2 + 2;
                end
            end
        end
    end

    for x_idx = 0:Lx-1
        for y_idx = 0:Ly -2
            % vertical bond
            nn_bond_idx_map(count1) = count2;
            count1 = count1 + 1;
            count2 = count2 + 1;
        end
        count2 = count2 + Ly - 2;  % skip NNN J2 interaction
    end
end


% ==== Plot Bond Energy ==== %
for i = 1:nn_bond_num
    x1_coor = nn_bond_coor_info(i, 1);
    y1_coor = nn_bond_coor_info(i, 2);
    x2_coor = nn_bond_coor_info(i, 3);
    y2_coor = nn_bond_coor_info(i, 4);

    bond_energy = bond_energys(nn_bond_idx_map(i));
    if(bond_energy > 0)
        bond_color = postive_bond_color;
    else
        bond_color = minus_bond_color;
    end
    line([x1_coor,x2_coor],[y1_coor,y2_coor],'color',bond_color,'linewidth',abs(bond_with_scale* bond_energy));  hold on;

    if(mark_energy_value)
        T = text((x1_coor+x2_coor)/2, (y1_coor + y2_coor)/2, [num2str(bond_energys(i))],'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

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
