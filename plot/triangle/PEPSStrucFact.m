Ly = 6;
Lx = 15;
J2 = 0;
Dpeps = 8;
Db = 8;
N = Lx * Ly;
auto_correlation_data_len=20;
site_num = Ly * Lx ;
num_points = floor(Lx / 2);
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

corr_data = fread(file_id, num_points * 6, 'double');
struc_factor_data = fread(file_id, N * N, 'double');
struc_factor_err = fread(file_id, N * N, 'double');
fclose(file_id);



