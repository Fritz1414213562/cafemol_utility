// This class is for the analysis on some dcd file.

//-------------------------------------------------------------------------------
// include guard
#pragma once
//-------------------------------------------------------------------------------
// includes

#include"DCDParser.hpp"
#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<array>
#include<numeric>
#include<algorithm>
#include<iterator>

//-------------------------------------------------------------------------------
// namespace
namespace caf {

//-------------------------------------------------------------------------------
// header

class DCDAnalyzer : public DCDParser {

private :
	std::array<std::array<float, 3>, 3> RodriguesRot(const std::array<float, 3>& norm, float ang);

    int frame_num;
    int atom_num;

    void read_num_frame_and_atom();

	std::array<float, 3> convertVec2Unit(const std::array<float, 3>& vec);

	std::array<float, 3> calcRelativeVec(const std::array<std::vector<float>, 3>& xyzs, int atom_id_begin, int atom_id_end);

	int searchNativeContact(const std::array<std::vector<float>, 3>& xyzs, const std::array<int, 2>& search_chain, const std::array<int, 2>& target_chain, const float cutoff);

	int searchNativeContact(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& search_chain, const std::vector<int>& target_chain, const float cutoff);


public :
    DCDAnalyzer();
    DCDAnalyzer(const std::string input_file_name);

	template<typename T, std::size_t N>
	std::array<T, N> calcDot(const std::array<std::array<T, N>, N>& lhs, const std::array<T, N>& rhs) {
		std::array<T, N> result;
		for (std::size_t idim = 0; idim < N; ++idim) {
			std::array<T, N> inner_pro;
			inner_pro = lhs[idim] * rhs;
			result[idim] = std::accumulate(inner_pro.begin(), inner_pro.end(), 0.0);
		}

		return result;
	}

	std::array<float, 3> Outer_pro(const std::array<float, 3>& vec1, const std::array<float, 3>& vec2);
    std::array<std::vector<std::array<float, 3>>, 3> makeVector(const std::array<int, 3>& atom_id_begins, const std::array<int, 3>& atom_id_ends);
	std::vector<std::array<float, 3>> makeVector(int atom_id_begin, int atom_id_end);
	std::array<std::array<float, 3>, 3> adjustRotMatrix(const std::array<float, 3>& reference_vec, const std::array<float, 3>& current_vec);

	std::vector<int> searchNativeContact(const std::array<int, 2>& search_chain, const std::array<int, 2>& target_chain, const float cutoff);

	std::vector<int> searchNativeContact(const std::array<int, 2>& search_chain, const std::array<int, 2>& target_chain, std::vector<int>& except_indices, const float cutoff);
};

// inline function---------------------------------------------------------------------------------


inline std::array<float, 3> caf::DCDAnalyzer::Outer_pro(const std::array<float, 3>& vec1, const std::array<float, 3>& vec2) {

	std::array<float, 3> result;
	result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

	return result;
}


inline std::array<std::array<float, 3>, 3> caf::DCDAnalyzer::RodriguesRot(const std::array<float, 3>& norm, float ang) {

	std::array<std::array<float, 3>, 3> result;
	for (std::size_t idim = 0; idim < 3; ++idim) {
		result[idim][idim] = std::cos(ang) + norm[idim] * norm[idim] * (1 - std::cos(ang));
	}
	result[0][1] = norm[0] * norm[1] * (1 - std::cos(ang)) - norm[2] * std::sin(ang);
	result[0][2] = norm[2] * norm[0] * (1 - std::cos(ang)) + norm[1] * std::sin(ang);
	result[1][0] = norm[0] * norm[1] * (1 - std::cos(ang)) + norm[2] * std::sin(ang);
	result[1][2] = norm[1] * norm[2] * (1 - std::cos(ang)) - norm[0] * std::sin(ang);
	result[2][0] = norm[2] * norm[0] * (1 - std::cos(ang)) - norm[1] * std::sin(ang);
	result[2][1] = norm[2] * norm[1] * (1 - std::cos(ang)) + norm[0] * std::sin(ang);

	return result;
}

}
