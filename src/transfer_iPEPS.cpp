//
// Created by haoxinwang on 03/11/2023.
//  No symmetry tensor

#include <iostream>
#include <fstream>
#include <vector>

#include "params_parser.h"
#include "gqdouble.h"
#include "gqpeps/two_dim_tn/peps/square_lattice_peps.h" //TPS, SquareLatticePEPS
#include "gqpeps/two_dim_tn/tps/split_index_tps.h"


using namespace gqten;
using namespace gqpeps;

void LoadTenData(const std::string &filename,
                 const size_t data_size,
                 double *data);

int main(int argc, char *argv[]) {
  VMCUpdateParams params(argv[1]);
  const size_t Lx = params.Lx;
  const size_t Ly = params.Ly;
  const size_t N = 3 * Lx * params.Ly;
  size_t Db = params.Db_max;
  const IndexT pb_out = IndexT({
                                   QNSctT(U1QN(0), 2)},
                               gqten::GQTenIndexDirType::OUT
  );
  const IndexT pb_in = gqten::InverseIndex(pb_out);

  const IndexT vb_out = IndexT({
                                   QNSctT(U1QN(0), Db)},
                               gqten::GQTenIndexDirType::OUT
  );
  const IndexT vb_in = gqten::InverseIndex(vb_out);
  const IndexT trvb_out = IndexT({
                                     QNSctT(U1QN(0), 1)},
                                 gqten::GQTenIndexDirType::OUT
  );
  const IndexT trvb_in = gqten::InverseIndex(trvb_out);

  Tensor Ba({pb_out, vb_out, vb_in});
  Ba({0, 0, 0}) = 1.0;
  Tensor Bb({pb_out, vb_in, vb_out});
  Bb({0, 0, 0}) = 1.0;
  Tensor Td({vb_out, vb_out, vb_in});
  Td({0, 0, 0}) = 1.0;
  Tensor Bc({pb_out, vb_in, vb_out});
  Bc({0, 0, 0}) = 1.0;
  Tensor Tu({vb_in, vb_in, vb_out});
  Tu({0, 0, 0}) = 1.0;
  Tensor filter_in({trvb_in, vb_out});
  filter_in({0, 0}) = 1.0;
  Tensor filter_out({trvb_out, vb_in});
  filter_out({0, 0}) = 1.0;

  std::string filename;
  std::string tensor_path = "../ipeps_tensors/D" + std::to_string(Db) + "/";

  filename = tensor_path + "IPESS_Ba_D" + std::to_string(Db) + ".dat";
  LoadTenData(filename, Ba.GetActualDataSize(), Ba.GetRawDataPtr());
  filename = tensor_path + "IPESS_Bb_D" + std::to_string(Db) + ".dat";
  LoadTenData(filename, Bb.GetActualDataSize(), Bb.GetRawDataPtr());
  filename = tensor_path + "IPESS_Bc_D" + std::to_string(Db) + ".dat";
  LoadTenData(filename, Bc.GetActualDataSize(), Bc.GetRawDataPtr());
  filename = tensor_path + "IPESS_Td_D" + std::to_string(Db) + ".dat";
  LoadTenData(filename, Td.GetActualDataSize(), Td.GetRawDataPtr());
  filename = tensor_path + "IPESS_Tu_D" + std::to_string(Db) + ".dat";
  LoadTenData(filename, Tu.GetActualDataSize(), Tu.GetRawDataPtr());

  Tensor tmp0, tmp1, tmp2, square_ten;
  Contract(&Tu, {0}, &Bc, {2}, &tmp0);
  Contract(&tmp0, {3}, &Td, {0}, &tmp1);
  Contract(&tmp1, {3}, &Bb, {1}, &tmp2);
  Contract(&tmp2, {3}, &Ba, {1}, &square_ten);
  square_ten.Transpose({0, 1, 4, 6, 2, 3, 5});
  /**
 * square_ten up to now
 *           3
 *           |
 * 0--combined_tps_ten--2
 *           |
 *           1
 *  and physical indexes:
 *    4: left_upper site
 *    5: lower site
 *    6: right site
 */

  size_t combined_phy_dim = 8;
  SplitIndexTPS<TenElemT, U1QN> split_idx_tps(Ly, Lx, combined_phy_dim);
  size_t phy_dim = 2;
  for (size_t row = 0; row < Ly; row++) {
    for (size_t col = 0; col < Lx; col++) {
      Tensor local_ten = square_ten;
      if (row == 0) {
        Tensor tmp;
        Contract(&local_ten, {3}, &filter_in, {1}, &tmp);
        local_ten = tmp;
        local_ten.Transpose({0, 1, 2, 6, 3, 4, 5});
      } else if (row == Ly - 1) {
        Tensor tmp;
        Contract(&local_ten, {1}, &filter_out, {1}, &tmp);
        local_ten = tmp;
        local_ten.Transpose({0, 6, 1, 2, 3, 4, 5});
      }

      if (col == 0) {
        Tensor tmp;
        Contract(&filter_in, {1}, &local_ten, {0}, &tmp);
        local_ten = tmp;
      } else if (col == Lx - 1) {
        Tensor tmp;
        Contract(&local_ten, {2}, &filter_out, {1}, &tmp);
        local_ten = tmp;
        local_ten.Transpose({0, 1, 6, 2, 3, 4, 5});
      }

      Tensor combined_tps_ten = local_ten; //use consistent name with the solver file.
      using QNT = U1QN;
      std::vector<Index<QNT>> indexes = combined_tps_ten.GetIndexes();
      std::vector<Index<QNT>> phy_indexes_out(indexes.begin() + 4, indexes.end());
      std::vector<Index<QNT>> phy_indexes_in(3);
      for (size_t i = 0; i < 3; i++) {
        phy_indexes_in[i] = InverseIndex(phy_indexes_out[i]);
      }
      split_idx_tps({row, col}) = std::vector<Tensor>(combined_phy_dim);
      for (size_t dim = 0; dim < combined_phy_dim; dim++) {
        Tensor proj_ten(phy_indexes_in);
        proj_ten({dim % phy_dim, (dim / phy_dim) % (phy_dim), (dim / phy_dim / phy_dim) % (phy_dim)}) = 1.0;
        Contract(&combined_tps_ten, {4, 5, 6}, &proj_ten, {0, 1, 2}, &split_idx_tps({row, col})[dim]);
        if (split_idx_tps({row, col})[dim].GetQNBlkNum() == 0) {
          std::cout << "warning: Site (" << row << "," << col << "), dim = " << dim << " projected tensor is empty."
                    << std::endl;
        }
      }
    }
  }

  split_idx_tps.Dump();
  return 0;
}

void LoadTenData(const std::string &filename,
                 const size_t data_size,
                 double *data) {// data is return, should alloc memory outside the function
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cout << "Failed to open file." << std::endl;
    exit(1);
  }

  int index = 0;
  while (file.eof() == false && index < data_size) {
    file >> data[index];
    index++;
  }

//  if (index != data_size) {
//    std::cout << "data_size : " << data_size << ", file data size : " << index << std::endl;
//    exit(1);
//  }
  file.close();

//  for(int i = 0; i < index; i++) {
//    std::cout << data[i] << std::endl;
//  }
  return;
}