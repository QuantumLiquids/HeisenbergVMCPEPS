/*
 * File Name: gqdouble.h
 * Description: Define the types and constant in Heisenberg PEPS project
 * Created by Hao-Xin on 2021/10/15.
 *
 */

#ifndef HUBBARD_SRC_GQDOUBLE_H
#define HUBBARD_SRC_GQDOUBLE_H

#include "gqten/gqten.h"


using TenElemT = gqten::GQTEN_Double;
using gqten::special_qn::U1QN;
using gqten::GQTensor;

using QNSctT = gqten::QNSector<U1QN>;
using IndexT = gqten::Index<U1QN>;
using Tensor = GQTensor<TenElemT, U1QN>;
const U1QN qn0 = U1QN(0); //N(particle number), Sz
#ifdef U1SYM
const IndexT pb_out = IndexT({QNSctT(U1QN(1), 1), QNSctT(U1QN(-1), 1)},
                             gqten::GQTenIndexDirType::OUT
);
#else
const IndexT pb_out = IndexT({
                                 QNSctT(U1QN({QNCard("Sz", U1QNVal(0))}), 2)},
                             gqten::GQTenIndexDirType::OUT
);
#endif
const IndexT pb_in = gqten::InverseIndex(pb_out);

#endif //HUBBARD_SRC_GQDOUBLE_H
