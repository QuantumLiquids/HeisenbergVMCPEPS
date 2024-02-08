Lx = 6;
Ly = 24;
J2 = 0;
Dpeps = 10;
Db = 10;

auto_correlation_data_len=20;
site_num = Ly * Lx ;
if(J2 == 0)
    filename1 = ['../../data/triangle_one_point_functions', num2str(Ly),'x', num2str(Lx),'D', num2str(Dpeps),'-',num2str(Db)];
    filename2 = ['../../data/triangle_one_point_functions', num2str(Ly),'x', num2str(Lx), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)];
    if(exist(filename2,"file"))
        file_id = fopen(filename2,'rb');
    else 
        file_id = fopen(filename1,'rb');
    end
else
    file_id = fopen(['../../data/triangle_one_point_functions', num2str(Ly),'x', num2str(Lx), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
end
sz = fread(file_id, site_num, 'double');
sz_err = fread(file_id, site_num, 'double');
sz_auto_corr = fread(file_id, auto_correlation_data_len, 'double');
fclose(file_id);


if ~isempty(sz_auto_corr)
    plot(1:auto_correlation_data_len, sz_auto_corr);
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

