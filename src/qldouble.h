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
#ifdef U1SYM
using QNT = qlten::special_qn::U1QN;
const QNT qn0 = QNT(0); //N(particle number), Sz
#else
using QNT = qlten::special_qn::TrivialRepQN;
const QNT qn0 = QNT(); //N(particle number), Sz
#endif
using qlten::QLTensor;

using QNSctT = qlten::QNSector<QNT>;
using IndexT = qlten::Index<QNT>;
using Tensor = QLTensor<TenElemT, QNT>;

#ifdef U1SYM
const IndexT pb_out = IndexT({QNSctT(QNT(1), 1), QNSctT(QNT(-1), 1)},
                             qlten::TenIndexDirType::OUT
);
#else
const IndexT pb_out = IndexT({QNSctT(qn0, 2)},
                             qlten::TenIndexDirType::OUT
);
#endif
const IndexT pb_in = qlten::InverseIndex(pb_out);

std::string peps_path = "peps";

#endif //SPIN_QLDOUBLE_H
