#include <fstream>

int main() {  
  std::ofstream outfile;

  outfile.open("test.txt", std::ios_base::app); // append instead of overwrite
  outfile << "Data"; 
  return 0;
}