clear;
Ly = 10;
Lx = 28;
J2 = 0.5;
Dpeps = 8;
Db = 24;


auto_correlation_data_len=20;
bond_num =  Lx * (Ly-1)+ (Lx-1) * Ly ;
site_num = Ly * Lx ;
if(J2 == 0)
    file_id = fopen(['../../data/square_one_point_functions', num2str(Lx),'x', num2str(Ly),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
else
    file_id = fopen(['../../data/square_one_point_functions', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
end
sz = fread(file_id, site_num, 'double');
sz_err = fread(file_id, site_num, 'double');
sz_auto_corr = fread(file_id, auto_correlation_data_len, 'double');
fclose(file_id);


if ~isempty(sz_auto_corr)
    semilogy(1:auto_correlation_data_len, sz_auto_corr);
    set(gca,'fontsize',24);
    set(gca,'linewidth',1.5);
    set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
    xlabel('$\Delta t$','Interpreter','latex');
    ylabel('Spin Configuration Auto-correlation','Interpreter','latex');
    set(get(gca,'XLabel'),'FontSize',24);
    set(get(gca,'YLabel'),'FontSize',24);
else
    error("No data found for spin auto-correlation!");
end

