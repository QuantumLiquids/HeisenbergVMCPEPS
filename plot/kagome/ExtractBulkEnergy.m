function [inner_bond_energy]=ExtractBulkEnergy(bond_energy_data, Ly, Lx, deep)
if nargin < 4
    deep = 1;
end

total_bond_num =  Lx*Ly*3 + Lx *(Ly-1)+ (Lx-1)*Ly + (Lx-1)*(Ly-1);
tri_bond_num = Lx*Ly*3;
vertical_bond_num = Lx *(Ly-1);
horizontal_bond_num = (Lx-1)*Ly;
diagonal_bond_num = (Lx-1)*(Ly-1);

vertical_bond_energy = bond_energy_data(end - vertical_bond_num+1:end);
vertical_bond_energy2=reshape(vertical_bond_energy, Lx,[])';
inner_vertical_bond_energy = vertical_bond_energy2(2:end-1, 2:end-1);
remain_bond_energy = bond_energy_data(1: end - vertical_bond_num);

data_segment = Lx * 3 + (Lx-1)*2;

remain_bond_energy2=remain_bond_energy(1: end - (3*Lx+(Lx-1)));
remain_bond_energy3=reshape(remain_bond_energy2, data_segment,[])';

A=remain_bond_energy3(1:end, 1:3*Lx+(Lx-1));
B=remain_bond_energy3(1:end, 3*Lx+Lx:end );
inner_bond_energy2=A(2:end-1, 4:end-4);
inner_bond_energy3=B(1:end-1, 1:end-1);

inner_bond_energy=[inner_bond_energy2(:); inner_bond_energy3(:); inner_vertical_bond_energy(:)];
end