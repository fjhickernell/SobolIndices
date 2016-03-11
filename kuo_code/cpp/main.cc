#include <iostream>
#include <fstream>
#include <vector>
#include "rref.h"


int main(){
  ifstream infile("joe-kuo-old.1111",ios::in);
  if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
    exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');

  for (unsigned j=1;j<=D-1;j++) {

    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
}




}


/*
    bool M1[3][3] = {   {1, 0, 1},
                {0, 1, 1},
                {1, 1, 0}   };
    BooleanMatrix booleanMatrix1(M1, 3, 3);
    std::cout << booleanMatrix1.getRank() << std::endl;

    bool M2[4][4] = {   {1,1,1,0},
                {0,1,1,0},
                {0,1,0,0},
                {1,1,1,1}   };
    BooleanMatrix booleanMatrix2(M2, 4, 4);
    std::cout << booleanMatrix2.getRank() << std::endl;

    std::vector< std::vector<bool> > MM(5);
    for (unsigned i = 0; i < 5; i++){
        for (unsigned j = 0; j < 5; j++){
                MM[i].push_back((i*j)%2);
        }
    }
    for (unsigned i = 0; i < 5; i++){
        for (unsigned j = 0; j < 5; j++){
                std::cout << MM[i][j] << " ";
        }
    std::cout << std::endl;
    }
    MM[4][4] = 1;
    BooleanMatrix booleanMatrix3(MM, 5, 5);
    std::cout << booleanMatrix3.getRank() << std::endl;
*/

