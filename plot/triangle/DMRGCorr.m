Ly = 8;
Lx = 16;
N = Lx * Ly;
J2 = 0.125;
D = 15000;
mark_energy_value = 0;
energy_text_size = 10;

bond_with_scale = 10;
marker_color =  [255,158,002] / 256; % Marker color for both positive and negative values
% marker_color =  [019, 103, 131]/256;
marker_face_color = marker_color; % Filled marker color

corr_file1 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx), 'J2', num2str(J2), 'D', num2str(D), 'corrzz.json'];
corr_data1 = jsondecode(fileread(corr_file1));
corr_file2 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx), 'J2', num2str(J2), 'D', num2str(D), 'corrpm.json'];
corr_data2 = jsondecode(fileread(corr_file2));
corr_file3 = ['../../data/triangle_dmrg', num2str(Ly), 'x', num2str(Lx), 'J2', num2str(J2), 'D', num2str(D), 'corrmp.json'];
corr_data3 = jsondecode(fileread(corr_file3));

DeltaX = zeros(1, numel(corr_data1));
SpinCorr = zeros(1, numel(corr_data1));
for i = 1:numel(corr_data1)
    site1_idx = corr_data1{i}{1}(1);
    site2_idx = corr_data1{i}{1}(2);
    SpinCorr(i) = corr_data1{i}{2} + 1/2 * corr_data2{i}{2} + 1/2 * corr_data3{i}{2};

    [x1_idx, y1_idx, x1_coor, y1_coor] = DMRGIdx2Coor(Ly, site1_idx);
    [x2_idx, y2_idx, x2_coor, y2_coor] = DMRGIdx2Coor(Ly, site2_idx);
    DeltaX(i) = x2_coor - x1_coor;
end

positiveIndices = SpinCorr >= 0;
negativeIndices = SpinCorr < 0;

h1=loglog(DeltaX(positiveIndices), abs(SpinCorr(positiveIndices)), 'o','Color', marker_color); hold on;
loglog(DeltaX(negativeIndices), abs(SpinCorr(negativeIndices)), 'o', 'Color', marker_color, 'MarkerFaceColor',  marker_face_color);


set(gca, 'fontsize', 24);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2); % Set line width 1.5 pounds
xlabel('$r$', 'Interpreter', 'latex');
ylabel('$\langle \bf S_i \cdot \bf S_j\rangle$', 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'FontSize', 24);
set(get(gca, 'YLabel'), 'FontSize', 24);