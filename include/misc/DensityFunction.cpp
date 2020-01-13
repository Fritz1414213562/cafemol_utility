#include"DensityFunction.hpp"


void cafemol::DensityFunction::set_DataRange(const float& data_left_lim, const float& data_right_lim){
	data_range[0] = data_left_lim;
	data_range[1] = data_right_lim;
}


std::array<std::vector<float>, 2> cafemol::DensityFunction::make_Density(const std::vector<float>& vec, const float bin_range) {

	std::array<std::vector<float>, 2> result;
	
	float data_left_lim;
	if (data_range.empty()) {
		std::vector<float>::const_iterator vec_min_itr = std::min_element(vec.begin(), vec.end());
		data_left_lim = *vec_min_itr;
	}
	else {
		data_left_lim = data_range[0];
	}
	
	int count = 0;
	float current_point = data_left_lim;
	float vector_size = vec.size();

	for (const float& vec_value : vec) {
		if ((current_point <= vec_value) && ((current_point + bin_range) > vec_value)) ++count;

		else if (current_point + bin_range <= vec_value) {
			float density = static_cast<float>(count) / vector_size;
			result[0].push_back(current_point);
			result[1].push_back(density);
			current_point += bin_range;
			count = 1;
		}
	}
	float density = static_cast<float>(count) / vector_size;
	result[0].push_back(current_point);
	result[1].push_back(density);

	return result;
}
