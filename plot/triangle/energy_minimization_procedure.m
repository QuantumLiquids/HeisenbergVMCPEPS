Lx = 8;
Ly = 8;
J2 = 0.11;
Dpeps = 10;


datadir = '../../data/';
filename = [datadir, 'energy_decreasing_triangle_',num2str(Ly),'x', num2str(Lx),'J20_11D', num2str(Dpeps), '.dat'];
energy_values = importdata(filename);

en_dmrg = -31.3887777147400335; % 8x8, D=12000
% 6x24 D=12000, en_dmrg = -71.1587492695341979;

h0=plot(energy_values, 'o','LineWidth', 2); hold on;
h1=yline(en_dmrg, 'b--');

l=legend([h0, h1], {'PEPS, $D=10$',  ...
    'DMRG $D=12000$'});
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',18);
set(l,'Location','NorthEast');


xlabel('Iteration','Interpreter','latex')
ylabel('total energy','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
set(gcf, 'position',  [100, 100, 1300, 600]);

hold off;
inset_axes = axes('Position', [0.3, 0.5, 0.35, 0.35]);

plot((energy_values - en_dmrg) / abs(en_dmrg), 'o', 'LineWidth', 2);hold on;

set(gca, 'YScale', 'log');
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2);
ylabel('energy relative error','Interpreter','latex')
box on;



figure_directory = '../../figure/triangle';
figure_name_eps = 'EnOpt8x8J2-0_11D10.eps';
figure_path = fullfile(figure_directory, figure_name_eps);
saveas(gcf, figure_path, 'epsc');
disp(['Energy deceasing of the heisenberg model: ', figure_path]);