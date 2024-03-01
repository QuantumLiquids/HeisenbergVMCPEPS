Ly = 8;
Lx = 8;
N = Lx * Ly;
J2 = 0;
D=12000;
mark_energy_value = 0;
energy_text_size=10;

bond_with_scale = 10;
postive_bond_color = [233, 196, 107]/256;
minus_bond_color = [042, 157, 142]/256;

corr_file1 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx),'J2', num2str(J2),'D', num2str(D),'zzsf.json'];
corr_data1 = jsondecode(fileread(corr_file1));
corr_file2 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx),'J2', num2str(J2),'D', num2str(D),'pmsf.json'];
corr_data2 = jsondecode(fileread(corr_file2));
corr_file3 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx),'J2', num2str(J2),'D', num2str(D),'mpsf.json'];
corr_data3 = jsondecode(fileread(corr_file3));

DeltaX = zeros(1, numel(corr_data1));
DeltaY = zeros(1, numel(corr_data1));
SpinCorr = zeros(1, numel(corr_data1));
for i = 1:numel(corr_data1)
    site1_idx = corr_data1{i}{1}(1);
    site2_idx = corr_data1{i}{1}(2);
    SpinCorr(i) = corr_data1{i}{2}  + 1/2 * corr_data2{i}{2} + 1/2 * corr_data3{i}{2};

    [x1_idx, y1_idx, x1_coor, y1_coor] = DMRGIdx2Coor(Ly,  site1_idx);
    [x2_idx, y2_idx, x2_coor, y2_coor] = DMRGIdx2Coor(Ly,  site2_idx);
    DeltaX(i) = x2_coor - x1_coor;
    DeltaY(i) = y2_coor - y1_coor;
end
% Calculate the spin structure factor in the Brillouin zone
ky_set = -2*pi/sqrt(3):0.02:2*pi/sqrt(3); % XC triangle lattice
kx_set = -4*pi/3:0.02:4*pi/3;
kx_num = numel(kx_set);
ky_num = numel(ky_set);
struct_factor = zeros(ky_num, kx_num);

for i = 1:numel(SpinCorr)
    struct_factor = struct_factor +  2 * SpinCorr(i) * cos( (kx_set *DeltaX(i) + ky_set' * DeltaY(i)));
end
struct_factor = (struct_factor + 3/4 * N)/ N/N;

% remove corner
for kx_idx = 1:kx_num
    kx = kx_set(kx_idx);
    for ky_idx = 1:ky_num
        ky = ky_set(ky_idx);
        if( abs(kx/(4*pi/3)) + abs(ky/(4*pi/sqrt(3))) > 1 )
            struct_factor(ky_idx, kx_idx) = NaN;
        end
    end
end

% Plotting the spin structure factor
figure;
imagesc(kx_set, ky_set,struct_factor);hold on;
colorbar;
T1=text(0,0,'$\Gamma$');
set(T1,'Interpreter','latex');set(T1,'Fontsize',32);
axis equal;
axis off;

set(gca,'fontsize',32);
set(gca,'linewidth',1.5);
%set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$M$','Interpreter','latex');
ylabel('$K$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',32);
set(get(gca,'YLabel'),'FontSize',32);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'linewidth',1.5);
% set(gcf,'position',[1000,1000, sqrt(3)*300, 2*300]);

% title('Spin Structure Factor');


