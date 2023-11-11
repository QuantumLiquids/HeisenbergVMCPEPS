file_id = fopen('statistic_summary','rb');
energy = fread(file_id, 1, 'double');

en_std = fread(file_id, 1, 'double');

sz = fread(file_id, 3 * 6*6, 'double');
plot(sz,'-o');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('site','Interpreter','latex');
ylabel('$S_z$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);