#include "gqmps2/gqmps2.h"
#include "myutil.h"
using gqmps2::kMpsPath;
using gqmps2::kMpsTenBaseName;
using gqmps2::kGQTenFileSuffix;

//number of mps file in default mps path("./mps")
size_t GetNumofMps(){
    size_t NumberOfMpsFile = 0;
    for(NumberOfMpsFile=0; NumberOfMpsFile < 1e5; NumberOfMpsFile++ ){
        std::string file;
        file = kMpsPath + "/" + kMpsTenBaseName + std::to_string(NumberOfMpsFile) + "." + kGQTenFileSuffix;
        std::ifstream ifs(file, std::ifstream::binary);
        if(ifs.good()){
            ifs.close();
        }else{
            break;
        }
    }
    return NumberOfMpsFile;
}


void Show(std::vector<size_t> v){
    for(auto iter = v.begin(); iter<v.end(); iter++){
        std::cout<< *iter<<",";
    }
    std::cout <<'\b' <<std::endl;
}


bool ParserBondDimension(int argc, char *argv[],
                         std::vector<size_t>& D_set) {
  int nOptionIndex = 1;
  std::string D_string;
  std::string arguement1 = "--D=";
  bool has_D_parameter(false);
  while (nOptionIndex < argc){
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0){
      D_string = &argv[nOptionIndex][arguement1.size()];
      has_D_parameter = true;
    }
    nOptionIndex++;
  }

  //split thread num list
  const char* split = ",";
  char *p;
  const size_t MAX_CHAR_LENTH = 1000;
  char D_char[MAX_CHAR_LENTH];
  for(size_t i =0;i< MAX_CHAR_LENTH;i++) {
    D_char[i] = 0;
  }

  strcpy(D_char, D_string.c_str() );

  p = strtok(D_char, split);
  while(p != nullptr){
    D_set.push_back( atoi(p) );
    p = strtok(nullptr, split);
  }

  return has_D_parameter;
}