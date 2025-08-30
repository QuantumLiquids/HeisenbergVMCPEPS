% Define the path to the data files
data_path = "two_point_function_samples";

% Get a list of all files in the data_path directory
file_list = dir(fullfile(data_path, 'sample*.csv'));

% Initialize a cell array to hold the data from each file
all_data = {};
combined_data = [];
% Loop through each file and read the data
for i = 1:length(file_list)
    % Get the full path of the current file
    file_name = fullfile(data_path, file_list(i).name);

    % Read the CSV file
    data = csvread(file_name);

    % Store the data in the cell array
    all_data{end+1} = data;
    combined_data = [combined_data; data];
end

% 8x8 
distance_max = 4;
distance_to_plot = 4;
corr_data = combined_data(:,distance_to_plot) ... \\Sz*Sz
            + 0.5 * sum(combined_data(:,[distance_to_plot + distance_max, distance_to_plot + 2 * distance_max]), 2);

% Plot a histogram of the column vector
figure;
% corr_data = corr_data(corr_data>0);
histogram(corr_data+0.25);%, 'Normalization', 'pdf'); % Normalize to show probability density

set(gca,'fontsize',18);
set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('spin correlation(r=4) $\langle S_i\cdot S_j\rangle$','Interpreter','latex');
ylabel('Sample Density');
yscale("log")
ylim([0.5, inf]);
xscale("log")
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'FontSize',18);

