Lx = 12;
Ly = 12;
Db = 8;

auto_correlation_data_len=20;
bond_num = Lx*Ly*3 + Lx *(Ly-1)+ (Lx-1)*Ly + (Lx-1)*(Ly-1);
file_id = fopen(['../data/kagome_statistic_summary', num2str(Ly),'x', num2str(Lx),'D', num2str(Db),'-iPEPS'],'rb');
energy = fread(file_id, 1, 'double');
en_std = fread(file_id, 1, 'double');
energy_auto_corr = fread(file_id, auto_correlation_data_len, 'double');
bond_energys = fread(file_id, bond_num, 'double');
sz = fread(file_id, 3* Lx*Ly, 'double');
spin_auto_corr=fread(file_id, auto_correlation_data_len, 'double');
% plot(sz,'-o');

fprintf(['energy : %.6f ', char(177),'%.3f\n'], energy, en_std);

e_site = 2 * mean(ExtractBulkEnergy(bond_energys, Ly, Lx));
fprintf(['energy persite(bulk) : %.6f.\n'], e_site);
fprintf(['energy persite(all) : %.6f.\n'], sum(bond_energys)/Lx/Ly/3);

x = spin_auto_corr - mean(sz.*sz);
semilogy(x,'-o'); hold on;
% plot(x(1:7),'-o'); hold on;

% symlog(gca,'y',-2)

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('site','Interpreter','latex');
ylabel('$A$','Interpreter','latex');
% ylabel('$S_z$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);