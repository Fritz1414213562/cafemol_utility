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
    std::array<int, 3> ids_begin = {308, 1312, 1638};
    //std::array<int, 3> ids_begin = {308, 1312, 1312};
    std::array<int, 3> ids_end = {311, 1802, 2128}; // chain C in above pdb file
    std::array<std::vector<std::array<float, 3>>, 3> vecs = dparser->makeVector(ids_begin, ids_end);

    // reference vector(unit vector)
    std::array<float, 3> reference_rel(vecs[0][0]);
	std::array<float, 3> reference_adj(vecs[1][0]);
	std::array<float, 3> reference_nor(vecs[2][0]);

    std::ofstream ofs;
    ofs.open(file_prefix + output_suffix);

	std::size_t ref_size = vecs[0].size();
	std::size_t vec_size = vecs[1].size();
	std::size_t nor_size = vecs[2].size();

	if ((ref_size != vec_size) || (ref_size != nor_size)) {
		std::cerr << "Error: The size of relative vec is not consistent to adjust vec or normal_vec" << std::endl;
		std::exit(1);
	}
    
	// rotate the vectors as the current adjust vector being consistent with reference one
    for (std::size_t iframe = 0; iframe < vec_size; ++iframe) {
		std::array<std::array<float, 3>, 3> Rot1 = dparser->adjustRotMatrix(reference_adj, vecs[1][iframe]);
		// rotate once
		std::array<float, 3> ivec = dparser->calcDot<float, 3>(Rot1, vecs[0][iframe]);
		std::array<float, 3> inor = dparser->calcDot<float, 3>(Rot1, vecs[2][iframe]);
		// rotate twice
		std::array<std::array<float, 3>, 3> Rot2 = dparser->adjustRotMatrix(reference_nor, inor);
		ivec = dparser->calcDot<float, 3>(Rot2, ivec);

        // calculate inner_product
		//std::array<float, 3> ivec2 = ivec * ivec;
		//float ivec_norm = std::sqrt(std::accumulate(ivec2.begin(), ivec2.end(), 0.0));
		//std::cout << ivec_norm << std::endl;
        std::array<float, 3> ref_and_ivec = reference_rel * ivec;
        const float inner_pro = std::accumulate(ref_and_ivec.begin(), ref_and_ivec.end(), 0.0);
        // calculate the angle between the reference vec and current vec
        float angle_inframe = 180 * std::acos(inner_pro) / PI;

        //ofs << inner_pro << std::endl;
        ofs << angle_inframe << std::endl;
    }
    ofs.close();

    return 0;
}

