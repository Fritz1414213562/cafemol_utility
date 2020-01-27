#include"DensityFunction.hpp"


//template<typename realT>
//void cafemol::DensityFunction<realT>::set_DataRange(const realT& data_left_lim, const realT& data_right_lim){
//	data_range[0] = data_left_lim;
//	data_range[1] = data_right_lim;
//}
//
//
//template<typename realT>
//std::array<std::vector<realT>, 2> cafemol::DensityFunction<realT>::make_Density(const std::vector<realT>& vec, const realT& bin_range) {
//
//	std::array<std::vector<realT>, 2> result;
//	
//	realT data_left_lim;
//	if (data_range.empty()) {
//		typename std::vector<realT>::const_iterator vec_min_itr = std::min_element(vec.begin(), vec.end());
//		data_left_lim = *vec_min_itr;
//	}
//	else {
//		data_left_lim = data_range[0];
//	}
//	
//	int count = 0;
//	realT current_point = data_left_lim;
//	realT vector_size = vec.size();
//
//	for (const realT& vec_value : vec) {
//		if ((current_point <= vec_value) && ((current_point + bin_range) > vec_value)) ++count;
//
//		else if (current_point + bin_range <= vec_value) {
//			realT density = static_cast<realT>(count) / vector_size;
//			result[0].push_back(current_point);
//			result[1].push_back(density);
//			current_point += bin_range;
//			count = 1;
//		}
//	}
//	realT density = static_cast<realT>(count) / vector_size;
//	result[0].push_back(current_point);
//	result[1].push_back(density);
//
//	return result;
//}
