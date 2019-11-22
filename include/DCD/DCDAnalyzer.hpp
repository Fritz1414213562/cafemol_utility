// This class is for the analysis on some dcd file.

//-------------------------------------------------------------------------------
// include guard
#ifndef DCD_ANALYZER_HPP
#define DCD_ANALYZER_HPP
//-------------------------------------------------------------------------------
// includes

#include"DCDParser.hpp"
#include"../misc/ContainerJudge.hpp"
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
namespace cafemol {

//-------------------------------------------------------------------------------
// header

class DCDAnalyzer : public DCDParser {

//-------------------------------------------------------------------------------
// public method

public :
    DCDAnalyzer();
    DCDAnalyzer(const std::string input_file_name);
	int get_AtomNumber();

	std::vector<int> getClosestResidue(const std::vector<int>& fixed_points_ids, const std::vector<int>& search_chain_ids);

	std::vector<int> getClosestResidue(const std::vector<int>& fixed_points_ids, const std::vector<int>& search_chain_ids, const float& cutoff);

	std::vector<std::array<int, 2>> getClosestResidue(const std::array<int, 2>& fixed_points_ids, const std::vector<int>& search_chain_ids, const float& cutoff);

	std::vector<int> get_AtomIDsNearPoint(const std::array<float, 3>& fixed_point_coordinate, const std::vector<int>& search_chain_ids, const float& cutoff);


	std::vector<std::array<int, 2>> get_1stNativeContactedResidue(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff);

	std::array<std::vector<int>, 2> get_NativeContactedRange(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff);

	std::vector<int> count_NativeContacts(const std::vector<std::array<int, 2>>& ini_contact_pairs, const float& cutoff);

	std::vector<int> count_NativeContacts(const std::vector<std::array<int, 2>>& ini_contact_pairs, const std::vector<float>& cutoff);

	std::vector<std::vector<int>> count_PDNSContacts(const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float, 3>>& parameters);

	std::vector<std::vector<std::array<int, 2>>> get_PDNSContactResidues(const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float, 3>>& parameters);

    std::vector<std::vector<int>> get_PDNSContactResidues(const std::vector<int>& pdns_pro_ids, const std::vector<int>& dna_ids, const float& cutoff);

	std::vector<int> get_NativeContacts(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff, const int& frame);

	std::vector<int> get_NativeContactsFromIni(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff);

    std::vector<std::vector<bool>> is_NativeContacts(const std::vector<int>& contact_search_ids, const std::vector<int>& contact_target_ids, const float& cutoff);


// ------------------------------------------------------------------------------
// private method
private :
	std::array<std::array<float, 3>, 3> RodriguesRot(const std::array<float, 3>& norm, float ang);

	std::array<float, 3> convertVec2Unit(const std::array<float, 3>& vec);

	std::array<float, 3> getVector(const std::array<std::vector<float>, 3>& xyzs, const int vec_id);

	std::array<float, 3> Outer_pro(const std::array<float, 3>& vec1, const std::array<float, 3>& vec2);

	float calcDistance(const std::array<float, 3>& vec1, const std::array<float, 3>& vec2) {

		std::array<float, 3> rel_vec = vec2 - vec1;
		std::array<float, 3> rel_vec2 = rel_vec * rel_vec;
		float result = std::sqrt(std::accumulate(rel_vec2.begin(), rel_vec2.end(), 0.0));

		return result;
	}

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

	template<typename Type, std::size_t N>
	Type calcDot(const std::array<Type, N>& lhs, const std::array<Type, N>& rhs) {
		Type result = 0.0;
		for (std::size_t idim = 0; idim < N; ++idim) {
			result += lhs[idim] * rhs[idim];
		}
		return result;
	}

	template<typename Type, std::size_t N>
	Type calc_Norm(const std::array<Type, N>& vec) {
		Type result = 0.0;
		for (std::size_t idim = 0; idim < N; ++idim) {
			result += vec[idim] * vec[idim];
		}
		result = std::sqrt(result);
		return result;
	}

	std::array<float, 3> calcCenter(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& ids);


	int getClosestResidue(const std::array<std::vector<float>, 3>& xyzs, const std::array<float, 3>& fixed_points_vector, const std::vector<int>& search_chain_ids);

	int getClosestResidue(const std::array<std::vector<float>, 3>& xyzs, const std::array<float, 3>& fixed_points_vector, const std::vector<int>& search_chain_ids, const float& cutoff);

	bool is_NativeContact(const std::array<float, 3>& vec1, const std::array<float, 3>& vec2, const float& cutoff);

	std::array<int, 2> get_1stNativeContactedResidue(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff);

	std::array<int, 2> get_lastNativeContactedResidue(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff);

	int count_NativeContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& ini_contact_pairs, const float& cutoff);

    int count_NativeContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& ini_contact_pairs, const std::vector<float>& cutoff);

	std::vector<int> count_PDNSContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float, 3>>& parameters);

	std::vector<std::array<int, 2>> get_PDNSContactResidues(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float, 3>>& parameters);

    std::vector<bool> is_NativeContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& contact_search_ids, const std::vector<int>& contact_target_ids, const float& cutoff);

	const float cutoff_angle = (10.0 / 180.0) * acos(-1);
};


// inline function---------------------------------------------------------------------------------


inline std::array<float, 3> cafemol::DCDAnalyzer::Outer_pro(const std::array<float, 3>& vec1, const std::array<float, 3>& vec2) {

	std::array<float, 3> result;
	result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

	return result;
}


inline std::array<std::array<float, 3>, 3> cafemol::DCDAnalyzer::RodriguesRot(const std::array<float, 3>& norm, float ang) {

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

#endif /* DCD_ANALYZER_HPP */
