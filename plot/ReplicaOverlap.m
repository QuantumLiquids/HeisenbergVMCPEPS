Lx = 6;
Ly = 6;
Db = 8;

auto_correlation_data_len=20;
bond_num = Lx*Ly*3 + Lx *(Ly-1)+ (Lx-1)*Ly + (Lx-1)*(Ly-1);
thread_num = 56;
overlap_sum = [];
overlap_data_collect=[];
for thr = 0:thread_num -1
    filename = ['../data/kagome_replica_overlap', num2str(Ly), 'x',num2str(Lx),'D', num2str(Db),'/replica_overlap',num2str(thr)];

    overlap_data = load(filename);
    if(thr == 0)
        overlap_sum = overlap_data;
        overlap_data_collect = overlap_data;
    else
        overlap_sum = overlap_sum + overlap_data;
        overlap_data_collect = [overlap_data_collect, overlap_data];
    end
    % fprintf(['energy : %.6f ', char(177),'%.3f\n'], energy, en_std);
    % plot([1;overlap_data],'-x'); hold on;
end
overlap_ave = mean(overlap_data_collect,2);
semilogx(overlap_ave,'-o'); hold on;
% std_err = std(overlap_data_collect,0,2)
% semilogx(std_err,'-x');hold on;


yline(0, 'r--'); hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('MC sweep','Interpreter','latex');
ylabel('overlap','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);