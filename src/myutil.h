
#ifndef HONEYCOMBHEISENBERGJ1J2J3_SRC_MYUTIL_H
#define HONEYCOMBHEISENBERGJ1J2J3_SRC_MYUTIL_H
#include <cstdlib>
#include <vector>
#include <string>

size_t GetNumofMps();
void Show(std::vector<size_t> v);
bool ParserBondDimension(int argc, char *argv[],
                         std::vector<size_t> &D_set);
bool IsFileExist(const std::string &);

enum ExtensionDir {
  ROW,
  COL
};
#endif
