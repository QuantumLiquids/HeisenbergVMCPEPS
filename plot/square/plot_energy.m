Lx = 12;
Ly = 12;
J2 = 0;
Dpeps = 8;
Db = 8;

auto_correlation_data_len=20;
bond_num =  Lx * (Ly-1)+ (Lx-1) * Ly ;
site_num = Ly * Lx ;
if(J2 == 0)
    file_id = fopen(['../../data/square_energy_statistics', num2str(Ly),'x', num2str(Lx),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
else
    file_id = fopen(['../../data/square_energy_statistics', num2str(Ly),'x', num2str(Lx), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
end
energy = fread(file_id, 1, 'double');
en_std = fread(file_id, 1, 'double');
bond_energys = fread(file_id, bond_num, 'double');
energy_auto_corr = fread(file_id, auto_correlation_data_len, 'double');
fclose(file_id);

fprintf('Energy +- en_std: %f pm %f\n', energy, en_std);
% fprintf('Bond Energies : \n');

% fprintf(['energy persite(bulk) : %.6f.\n'], e_site);
% fprintf(['energy persite(all) : %.6f.\n'], sum(bond_energys)/Lx/Ly/3);

if ~isempty(energy_auto_corr)
    plot(1:auto_correlation_data_len, energy_auto_corr);
    set(gca,'fontsize',24);
    set(gca,'linewidth',1.5);
    set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
    xlabel('$\Delta t$','Interpreter','latex');
    ylabel('Energy Auto-correlation','Interpreter','latex');
    set(get(gca,'XLabel'),'FontSize',24);
    set(get(gca,'YLabel'),'FontSize',24);
else
    error("No data found for energy autocorrelation!");
end

