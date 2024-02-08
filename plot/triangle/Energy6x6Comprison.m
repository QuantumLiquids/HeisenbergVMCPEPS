% Read the data from the file
data = dlmread('energy6x6.txt', '\t', 1, 0);

% Extract the columns for J2, DMRGD6000, and PEPSD8
J2 = data(:, 1);
DMRGD6000 = data(:, 2);
PEPSD8 = data(:, 3);

% Plot the first figure
figure(1);
subplot(1, 2, 1);
plot(J2, DMRGD6000, '-o', 'LineWidth', 2);
hold on;
plot(J2, PEPSD8, '-x', 'LineWidth', 2);
hold off;

% Add labels and a legend to the first figure
l = legend('DMRG, $D=6000$', 'PEPS, $D=8$');
set(l,'Box','off');
set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','Best');
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2);
xlabel('$J_2$','Interpreter','latex');
ylabel('Total Energy','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

% Calculate the relative error between the two curves
relative_err = (PEPSD8 - DMRGD6000) ./ abs(DMRGD6000);

% Plot the second figure
subplot(1, 2, 2);
semilogy(J2, relative_err, '-s', 'LineWidth', 2);

% Add labels and a title to the second figure
xlabel('$J_2$', 'Interpreter', 'latex');

% Move the y-axis label to the right side
ylabel('Energy Relative Error', 'Interpreter', 'latex');

% Adjust the plot appearance
set(gca, 'fontsize', 24);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2);
set(get(gca, 'XLabel'), 'FontSize', 24);
set(get(gca, 'YLabel'), 'FontSize', 24);

% Adjust the spacing between the subplots
subplot(1, 2, 1);
pos = get(gca, 'Position');
pos(3) = pos(3) - 0.01;
set(gca, 'Position', pos);

subplot(1, 2, 2);
pos = get(gca, 'Position');
pos(1) = pos(1) + 0.01;
set(gca, 'Position', pos);