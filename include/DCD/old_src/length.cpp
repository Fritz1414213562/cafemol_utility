#include"DCDAnalyzer.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<memory>
#include<array>
#include<vector>
#include<numeric>
#include<cmath>
#include<typeinfo>

template<typename T, std::size_t N> 
inline std::array<T, N> operator*(const std::array<T, N> lhs, const std::array<T, N> rhs) {

    std::array<T, N> result;
    for (std::size_t iN = 0; iN < N; ++iN) {
        result[iN] = lhs[iN] * rhs[iN];
    }

    return result;
}


int main(int argc, char *argv[]) {

    // constant
    const float PI = 4 * atan(1);

    std::string file_name = argv[1];
	std::string output_name = argv[2];
	std::string file_prefix;
	std::string file_suffix;
	std::size_t file_index = 0;

	// judge the command line ---------------------------------------------------
	for (const char& str : file_name) {

		if (str == '.') {break;}

		else if (file_index < file_name.size()) {
			file_prefix += str;
			++file_index;
		}

		else {
			std::cerr << "Error: DCD file does not exist" << std::endl;
			std::exit(1);
		}
	}

	for (std::size_t idx = file_index; idx < file_name.size(); ++idx) {

		file_suffix += file_name[idx];
	}

	if (file_suffix != ".dcd") {
		std::cerr << "Error: This file does not have the suffix '.dcd'" << std::endl;
		std::exit(1);
	}

	else {
		std::cout << "The file '" << file_prefix << ".dcd' has read" << std::endl;
	}

    std::string input_suffix = ".dcd";
    std::string output_suffix = ".txt";
	// ---------------------------------------------------------------------------------
	// calculation

    std::unique_ptr<caf::DCDAnalyzer>dparser =  std::make_unique<caf::DCDAnalyzer>(file_prefix + input_suffix);
    int ids_begin = 1481;
    int ids_end = 995; 
    std::vector<std::array<float, 3>> vecs = dparser->makeVector(ids_begin, ids_end);

    // reference vector(unit vector)
    std::array<float, 3> reference_rel(vecs[0]);
	std::array<float, 3> ref2 = reference_rel * reference_rel;
	float ref_length = std::sqrt(std::accumulate(ref2.begin(), ref2.end(), 0.));

    std::ofstream ofs;
    ofs.open(output_name + output_suffix);

	// calculate the length
    for (const std::array<float, 3>& ivec : vecs) {

        std::array<float, 3> ivec2 = ivec * ivec;
        float ivec_length = std::sqrt(std::accumulate(ivec2.begin(), ivec2.end(), 0.0));
		float dlength = ivec_length - ref_length;

        //ofs << inner_pro << std::endl;
        ofs << dlength << std::endl;
    }
    ofs.close();

    return 0;
}

