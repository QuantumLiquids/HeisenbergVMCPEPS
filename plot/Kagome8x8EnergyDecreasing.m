e0D8=[-82.03496291, -81.89328610, -81.72985196, -81.87593464, -81.91994767, -81.94854546, -81.96776235, -81.98130512, -81.99068962, -81.99032571, -82.00979211, -81.99524703, -82.01573015, -82.01137906, -82.02208827, -82.02550178, -82.02487136, -82.03243577, -82.03488748, -82.03472566, -82.04393498, -82.02992805, -82.04254445, -82.04308069, -82.05339256, -82.03848862, -82.04802953, -82.05923634, -82.06211603, -82.04225070, -82.04947195, -82.05713468, -82.07071413, -82.07626396, -82.07328963, -82.07207495, -82.07912691, -82.08612838, -82.08646264, -82.08136615, -82.08767503, -82.08781033, -82.07049913, -82.09380423, -82.09666653, -82.08309535, -82.08952325, -82.08693912, -82.08884191, -82.09273296, -82.08809097, -82.09369396, -82.09368290, -82.09627742, -82.09276665, -82.09113315, -82.10576376, -82.10533722, -82.10870465, -82.10922403];
e0D10=[];

h0=plot(e0D8([1:end]), 'o','LineWidth', 2); hold on;

% h0=plot(e0D10([1,3:end]), 'o','LineWidth', 2); hold on;
% 
% T=text(200,0.03,['$4\times 4$, PEPS $D=8$', char(10),...
%                 'Boundary MPS $\chi=24$', char(10),...
%                 'Typical step length $\alpha = 0.2$', char(10),...
%                 'Typical Samples = 56 chains $\times 2000$.']);
% set(T,'Interpreter','latex');
% set(T,'Fontsize',18);

% 
% l=legend([h0, h1, h2,h3, h4,h5], {'SR, $D=8$', 'Gradient, $D=8$', ...
%     'Gradient, $D=8, \alpha = 0.01$', 'SR, $D=10, \alpha = 0.2$', ...
%     'SR, $D=10$', ... %$\alpha=0.3$ for initial and final stage, varies in middle stage
%     'SR, $D=12$'});
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',18);
% set(l,'Location','NorthEast');

xlabel('Iteration','Interpreter','latex')
% ylabel('Persite energy $e_0$','Interpreter','latex')
ylabel('energy relative error','Interpreter','latex')
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
set(gcf, 'position',  [100, 100, 1300, 600]);
hold off;

figure_directory = './';
figure_name_eps = 'EnergyDecreasingKagome6x6.eps';
figure_path = fullfile(figure_directory, figure_name_eps);
saveas(gcf, figure_path, 'epsc');
disp(['Energy deceasing of the heisenberg model: ', figure_path]);