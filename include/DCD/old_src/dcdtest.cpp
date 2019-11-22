#include"DCDParser.hpp"
#include<string>
#include<array>
#include<vector>
#include<memory>
#include<fstream>

int main() {

    std::string input_file_name = "./data/moveOut.dcd"
    std::unique_ptr<caf::DCDParser> parser = std::make_unique<caf::DCDParser>(input_file_name);
    std::array<std::vector<float>, 3> xyzs = parser->read_binary_file();
    std::string output_file_name = "./script/moveOut.txt"
    std::ofstream ofs(output_file_name);
    std::size_t atom_num = xyzs[0].size();
    for (std::size_t atom_i = 0; atom_i < atom_num; ++atom_i) {
        ofs << xyzs[0][atom_i] << " " << xyzs[1][atom_i] << " " << xyzs[2][atom_i] << std::end;
    }
    ofs.close()
    return 0;
}
