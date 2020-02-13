#ifndef SET_SIMILARITY_CALCULATION_HPP
#define SET_SIMILARITY_CALCULATION_HPP
#include<UtilFunc.hpp>
#include<array>
#include<vector>


namespace cafemol::library {


template<typename Size3_Vector>
const float calc_TotalElementSize(const std::vector<Size3_Vector>& vec) {
	float result = 0.0;
	for (const Size3_Vector& vec_elem : vec) {
		result += vec_elem[2];
	}
	return result;
}



template<typename Size3_Vector>
void overlie_Data(std::vector<Size3_Vector>& overlain_data, const std::vector<Size3_Vector>& input_data, const float& bin_width) {
	for (const Size3_Vector& input_datum : input_data) {
		bool is_overlain = false;
		std::size_t overlain_index = 0;
		for (Size3_Vector& overlain_datum : overlain_data) {
			if (is_on_same_pixel(input_datum[0], input_datum[1],
								 overlain_datum[0], overlain_datum[1], bin_width)) {
				is_overlain = true;
				break;
			}
			++overlain_index;
		}
		if (is_overlain) overlain_data[overlain_index][2] += input_datum[2];
		else overlain_data.push_back(input_datum);
	}
}






template<typename Size3_Vector>
const float calc_SimpsonSimilarity(const std::vector<Size3_Vector>& lhs, const std::vector<Size3_Vector>& rhs, const float& bin_width) {

	float elem_contrib_sum = 0.0;
	const float& lhs_total_elem_size = calc_TotalElementSize(lhs);
	const float& rhs_total_elem_size = calc_TotalElementSize(rhs);
	const float& smaller_total_elem_size = std::min(lhs_total_elem_size, rhs_total_elem_size);

	for (const Size3_Vector& lhs_elem : lhs) {
		bool is_same_elem = false;
		std::size_t same_elem_index = 0;
		for (const Size3_Vector& rhs_elem : rhs) {
			is_same_elem = is_on_same_pixel(lhs_elem[0], lhs_elem[1],
											rhs_elem[0], rhs_elem[1], bin_width);
			if (is_same_elem) break;
			++same_elem_index;
		}
		if (is_same_elem) elem_contrib_sum += std::min(lhs_elem[2], rhs[same_elem_index][2]);
	}

	return 1.0 - (elem_contrib_sum / smaller_total_elem_size);
}




// Dice Index


template<typename Size3_Vector>
const float calc_DiceSimilarity(const std::vector<Size3_Vector>& lhs, const std::vector<Size3_Vector>& rhs, const float& bin_width) {

	float elem_contrib_sum = 0.0;
	const float& total_contrib_size = calc_TotalElementSize(lhs) + calc_TotalElementSize(rhs);

	for (const Size3_Vector& lhs_elem : lhs) {
		bool is_same_elem = false;
		std::size_t same_elem_index = 0;
		for (const Size3_Vector& rhs_elem : rhs) {
			is_same_elem = is_on_same_pixel(lhs_elem[0], lhs_elem[1],
											rhs_elem[0], rhs_elem[1], bin_width);
			if (is_same_elem) break;
			++same_elem_index;
		}
		if (is_same_elem) elem_contrib_sum += std::min(lhs_elem[2], rhs[same_elem_index][2]);
	}

	return 1.0 - (2.0 * elem_contrib_sum / total_contrib_size);
}




// Jaccard Index


template<typename Size3_Vector>
const float calc_JaccardSimilarity(const std::vector<Size3_Vector>& lhs, const std::vector<Size3_Vector>& rhs, const float& bin_width) {

	float elem_contrib_sum = 0.0;
	const float& total_contrib_size = calc_TotalElementSize(lhs) + calc_TotalElementSize(rhs);

	for (const Size3_Vector& lhs_elem : lhs) {
		bool is_same_elem = false;
		std::size_t same_elem_index = 0;
		for (const Size3_Vector& rhs_elem : rhs) {
			is_same_elem = is_on_same_pixel(lhs_elem[0], lhs_elem[1],
											rhs_elem[0], rhs_elem[1], bin_width);
			if (is_same_elem) break;
			++same_elem_index;
		}
		if (is_same_elem) elem_contrib_sum += std::min(lhs_elem[2], rhs[same_elem_index][2]);
	}

	const float& union_contrib_size = total_contrib_size - elem_contrib_sum;

	return 1.0 - (elem_contrib_sum / union_contrib_size);
}








}



#endif /* SET_SIMILARITY_CALCULATION_HPP */
