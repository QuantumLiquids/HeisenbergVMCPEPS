clear;
Ly = 10;
Lx = 28;
J2 = 0.5;
Dpeps = 8;
Db = 24;

auto_correlation_data_len=20;
site_num = Ly * Lx ;
num_points = floor(Lx / 2);
if(J2 == 0)
    filename1 = ['../../data/square_two_point_functions', num2str(Lx),'x', num2str(Ly),'D', num2str(Dpeps),'-',num2str(Db)];
    filename2 = ['../../data/square_two_point_functions', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)];
    if(exist(filename2,"file"))
        file_id = fopen(filename2,'rb');
    else
        file_id = fopen(filename1,'rb');
    end
else
    file_id = fopen(['../../data/square_two_point_functions', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
end

data = fread(file_id, num_points * 6, 'double');
fclose(file_id);

sz_sz_corr = data(1:num_points);
s_plus_s_minus_corr = data(num_points+1:2*num_points);
s_minus_s_plus_corr = data(2*num_points+1:3*num_points);
sz_sz_corr_err = data(3*num_points+1:4*num_points);
s_plus_s_minus_corr_err = data(4*num_points+1:5*num_points);
s_minus_s_plus_corr_err = data(5*num_points+1:6*num_points);

% Plot the spin correlation
x = 1:num_points;  
% h1 = errorbar(x, sz_sz_corr, sz_sz_corr_err);
% hold on;
% h2 = errorbar(x, 0.5 * s_plus_s_minus_corr, 0.5 * s_plus_s_minus_corr_err);
% h3 = errorbar(x, 0.5 * s_minus_s_plus_corr, 0.5 * s_minus_s_plus_corr_err);
% hold off;

ss_corr = sz_sz_corr + 0.5 * s_plus_s_minus_corr + 0.5 * s_minus_s_plus_corr;
ss_corr_err = sqrt(sz_sz_corr_err.^2 + (0.5 * s_plus_s_minus_corr_err).^2 + (0.5 * s_minus_s_plus_corr_err).^2);

h = errorbar(x, abs(ss_corr), ss_corr_err, '-x'); hold on;
% errorbar(x, 1.5 * abs(sz_sz_corr), 1.5 * sz_sz_corr_err);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\Delta x$','Interpreter','latex');
ylabel('$|\langle \bf S_i \cdot \bf S_j\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);



