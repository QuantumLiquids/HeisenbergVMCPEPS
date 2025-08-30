Ly = 8;
Lx = 24;
J2 = 0.125;
Dpeps = 10;
Db = 20;

auto_correlation_data_len=20;
site_num = Ly * Lx ;
num_points = floor(Lx / 2);

marker_color =  [019, 103, 131]/256; % Marker color for both positive and negative values
marker_face_color = marker_color; % Filled marker color
marker_size = 8;
if(J2 == 0)
    filename1 = ['../../data/triangle_two_point_functions', num2str(Lx),'x', num2str(Ly),'D', num2str(Dpeps),'-',num2str(Db)];
    filename2 = ['../../data/triangle_two_point_functions', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)];
    if(exist(filename2,"file"))
        file_id = fopen(filename2,'rb');
    else 
        file_id = fopen(filename1,'rb');
    end
else
    file_id = fopen(['../../data/triangle_two_point_functions', num2str(Lx),'x', num2str(Ly), 'J2',num2str(J2),'D', num2str(Dpeps),'-',num2str(Db)],'rb');
end

corr_data = fread(file_id, num_points * 3, 'double');
struc_factor_data = fread(file_id, site_num * site_num, 'double');
corr_err_data = fread(file_id, num_points * 3, 'double');
struc_factor_err = fread(file_id, site_num * site_num, 'double');

fclose(file_id);

sz_sz_corr = corr_data(1:num_points);
s_plus_s_minus_corr = corr_data(num_points+1:2*num_points);
s_minus_s_plus_corr = corr_data(2*num_points+1:3*num_points);
sz_sz_corr_err = corr_err_data(1:num_points);
s_plus_s_minus_corr_err = corr_err_data(num_points+1:2*num_points);
s_minus_s_plus_corr_err = corr_err_data(2*num_points+1:3*num_points);

% Plot the spin correlation
x = 1:num_points;  
% h1 = errorbar(x, sz_sz_corr, sz_sz_corr_err);
% hold on;
% h2 = errorbar(x, 0.5 * s_plus_s_minus_corr, s_plus_s_minus_corr_err);
% h3 = errorbar(x, 0.5 * s_minus_s_plus_corr, s_minus_s_plus_corr_err);
% hold off;

ss_corr = sz_sz_corr + 0.5 * s_plus_s_minus_corr + 0.5 * s_minus_s_plus_corr;
ss_corr_err = sqrt(sz_sz_corr_err.^2 + (0.5 * s_plus_s_minus_corr_err).^2 + (0.5 * s_minus_s_plus_corr_err).^2);

positiveIndices = ss_corr >= 0;
negativeIndices = ss_corr < 0;


positive_h2 = errorbar(x(positiveIndices), abs(ss_corr(positiveIndices)), ss_corr_err(positiveIndices), 'diamond'); hold on;
negative_h = errorbar(x(negativeIndices), abs(ss_corr(negativeIndices)), ss_corr_err(negativeIndices), 'diamond');

plot(x, abs(ss_corr),'-.', 'Color',marker_color);

set(positive_h2, 'Color', marker_color);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
set(positive_h2, 'MarkerSize', marker_size);

set(negative_h, 'Color', marker_color);
set(negative_h, 'MarkerFaceColor',  marker_face_color);
set(negative_h, 'MarkerSize', marker_size);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('$\langle \bf S_i \cdot \bf S_j\rangle$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);



