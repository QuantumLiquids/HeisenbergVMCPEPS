
num_mc_chain = 100;
all_data =[];
for i = 1:num_mc_chain
    en_loc_data = load(['energy_raw_data\energy', num2str(i-1)]);
    all_data = [all_data; en_loc_data];
end
all_data = all_data(all_data<0);
h = histogram(all_data, 'NumBins', 100000);
ylim([min(0.1), max(h.Values)]);

set(gca, 'YScale', 'log');