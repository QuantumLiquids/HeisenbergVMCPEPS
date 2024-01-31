/*
 * File Name: qldouble.h
 * Description: Define the types and constant in Heisenberg PEPS project
 * Created by Hao-Xin on 2021/10/15.
 *
 */

#ifndef SPIN_QLDOUBLE_H
#define SPIN_QLDOUBLE_H

#include "qlten/qlten.h"

using TenElemT = qlten::QLTEN_Double;
using qlten::special_qn::U1QN;
using qlten::QLTensor;

using QNSctT = qlten::QNSector<U1QN>;
using IndexT = qlten::Index<U1QN>;
using Tensor = QLTensor<TenElemT, U1QN>;
const U1QN qn0 = U1QN(0); //N(particle number), Sz
#ifdef U1SYM
const IndexT pb_out = IndexT({QNSctT(U1QN(1), 1), QNSctT(U1QN(-1), 1)},
                             qlten::TenIndexDirType::OUT
);
#else
const IndexT pb_out = IndexT({
                                 QNSctT(U1QN(0), 2)},
                             qlten::TenIndexDirType::OUT
);
#endif
const IndexT pb_in = qlten::InverseIndex(pb_out);

std::string peps_path = "peps";

#endif //SPIN_QLDOUBLE_H
