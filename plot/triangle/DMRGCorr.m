Ly = 8;
Lx = 24;
J2 = 0.125;
D=15000;
mark_energy_value = 0;
energy_text_size=10;

bond_with_scale = 10;
postive_bond_color = [233, 196, 107]/256;
minus_bond_color = [042, 157, 142]/256;

corr_file1 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx),'J2', num2str(J2),'D', num2str(D),'corrzz.json'];
corr_data1 = jsondecode(fileread(corr_file1));
corr_file2 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx),'J2', num2str(J2),'D', num2str(D),'corrpm.json'];
corr_data2 = jsondecode(fileread(corr_file2));
corr_file3 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx),'J2', num2str(J2),'D', num2str(D),'corrmp.json'];
corr_data3 = jsondecode(fileread(corr_file3));

DeltaX = zeros(1, numel(corr_data1));
SpinCorr = zeros(1, numel(corr_data1));
for i = 1:numel(corr_data1)
    site1_idx = corr_data1{i}{1}(1);
    site2_idx = corr_data1{i}{1}(2);
    SpinCorr(i) = corr_data1{i}{2}  + 1/2 * corr_data2{i}{2} + 1/2 * corr_data3{i}{2};

    [x1_idx, y1_idx, x1_coor, y1_coor] = DMRGIdx2Coor(Ly,  site1_idx);
    [x2_idx, y2_idx, x2_coor, y2_coor] = DMRGIdx2Coor(Ly,  site2_idx);
    DeltaX(i) = x2_coor - x1_coor;
end
loglog(DeltaX, abs(SpinCorr) ,'-o'); hold on;

% axis tight;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\Delta x$','Interpreter','latex');
ylabel('$\langle \bf S_i \cdot \bf S_j\rangle$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);


