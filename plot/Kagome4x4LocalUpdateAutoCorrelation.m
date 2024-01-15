data = load('../data/kagome_energy_raw_data8x8D8-iPESS/energy2');
t_list = 1:100;
A_list = zeros(1,numel(t_list)); %auto correlation
mu = mean(data);
std_energ = std(data);
for i = 1:numel(t_list)
    t = t_list(i);
    A_list(i) = mean(data(1:end-t) .*data(1+t:end))-mu*mu;
end
plot(t_list, (A_list)/std_energ,'-o'); hold on;

fit_model = fit(t_list(1:5)', abs(A_list(1:5))', 'exp1');
% semilogy(t_list(1:5), fit_model(t_list(1:5)));
fprintf('Autocorrelation Length: %.2f\n', -1 / fit_model.b);

% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$t$','Interpreter','latex');
ylabel('$A(t)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
