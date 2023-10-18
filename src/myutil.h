
#ifndef HONEYCOMBHEISENBERGJ1J2J3_SRC_MYUTIL_H
#define HONEYCOMBHEISENBERGJ1J2J3_SRC_MYUTIL_H
#include <stdlib.h>
#include <vector>

size_t GetNumofMps();
void Show(std::vector<size_t> v);
bool ParserBondDimension(int argc, char *argv[],
                         std::vector<size_t>& D_set);

#endif
