% Define the path to the data files
data_path = "wave_function_amplitudes";

% Get a list of all files in the data_path directory
file_list = dir(fullfile(data_path, 'psi*'));

combined_data = [];
% Loop through each file and read the data
for i = 1:length(file_list)
    % Get the full path of the current file
    file_name = fullfile(data_path, file_list(i).name);

    % Read the CSV file
    data = load(file_name);

    combined_data = [combined_data; data];
end

histogram(abs(combined_data));%, 'Normalization', 'pdf'); % Normalize to show probability density

set(gca,'fontsize',18);
set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\Psi(S)$','Interpreter','latex');
ylabel('Sample Density');
yscale("log")
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'FontSize',18);

