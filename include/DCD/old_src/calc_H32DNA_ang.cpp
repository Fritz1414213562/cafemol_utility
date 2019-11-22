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


int main() {

    // constant
    const float PI = 4 * atan(1);

    std::string input_file_name = "test.dcd";
    std::string output_file_name = "test_angle.txt";

    std::unique_ptr<caf::DCDAnalyzer>dparser =  std::make_unique<caf::DCDAnalyzer>(input_file_name);
    const int base1_id = 308; // Chain A in 1kx5_601-200bp_CGmodel.pdb
    const int base2_id = 311; // chain C in 1kx5_601-200bp_CGmodel.pdb
    const int histone8mer_id_begin = 1199;
    const int histone8mer_id_end = 2178;
    std::vector<std::array<float, 3>> vecs = dparser->makeVector(base1_id, base2_id, histone8mer_id_begin, histone8mer_id_end);

    // reference vector(unit vector)
    std::array<float, 3> reference_vector(vecs[0]);

    std::ofstream ofs;
    ofs.open(output_file_name);
    
    for (const std::array<float, 3>& ivec : vecs) {

        // calculate inner_product
        std::array<float, 3> ref_and_ivec = reference_vector * ivec;
        const float inner_pro = std::accumulate(ref_and_ivec.begin(), ref_and_ivec.end(), 0.0);
        // calculate the angle between the reference vec and current vec
        float angle_inframe = 180 * std::acos(inner_pro) / PI;

        //ofs << inner_pro << std::endl;
        ofs << angle_inframe << std::endl;
    }
    ofs.close();

    return 0;
}

