Ly = 10;
Lx = 10;
J2 = 0;
Dpeps = 10;
Db = 10;
N = Lx * Ly;
auto_correlation_data_len=20;
site_num = Ly * Lx ;
num_points = floor(Lx / 2);
if(J2 == 0)
    filename1 = ['../../data/triangle_two_point_functions', num2str(Lx),'x', num2str(Ly),'D', num2str(Dpeps),'-',num2str(Db)];
    filename2 = ['../../data/triangle_two_point_functions', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)];
    if(exist(filename2,"file"))
        file_id = fopen(filename2,'rb');
    else 
        file_id = fopen(filename1,'rb');
    end
else
    file_id = fopen(['../../data/triangle_two_point_functions', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
end


corr_data = fread(file_id, num_points * 3, 'double');
struc_factor_data = fread(file_id, site_num * site_num, 'double');

corr_err_data = fread(file_id, num_points * 3, 'double');
struc_factor_err = fread(file_id, site_num * site_num, 'double');

fclose(file_id);

DeltaX = zeros(1, site_num * site_num);
DeltaY = zeros(1, site_num * site_num);
SpinCorr = 3 * struc_factor_data; % only sz
i = 1;
for site1_idx = 0:N-1
    [x1_idx, y1_idx, x1_coor, y1_coor] = PEPSIdx2Coor(Lx,  site1_idx);
    for site2_idx = 0:N-1
        [x2_idx, y2_idx, x2_coor, y2_coor] = PEPSIdx2Coor(Lx,  site2_idx);
        DeltaX(i) = x2_coor - x1_coor;
        DeltaY(i) = y2_coor - y1_coor;
        i = i+1;
    end
end


% Calculate the spin structure factor in the Brillouin zone
ky_set = -2*pi/sqrt(3):0.02:2*pi/sqrt(3); % XC triangle lattice
kx_set = -4*pi/3:0.02:4*pi/3;
kx_num = numel(kx_set);
ky_num = numel(ky_set);
struct_factor = zeros(ky_num, kx_num);

for i = 1:numel(SpinCorr)
    struct_factor = struct_factor + SpinCorr(i) * cos( (kx_set *DeltaX(i) + ky_set' * DeltaY(i)));
end
struct_factor = struct_factor/ N/N;

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
% T1=text(0,0,'$\Gamma$');
% set(T1,'Interpreter','latex');set(T1,'Fontsize',32);
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




