#ifndef DENSITY_FUNCTION_HPP
#define DENSITY_FUNCTION_HPP
#include<UtilContainer.hpp>
#include<ErrorMessage.hpp>
#include<cmath>
#include<vector>
#include<array>
#include<algorithm>
#include<tuple>

namespace cafemol {


template<typename InputType, typename YAxisType>
using HistogramData = std::tuple<InputType, YAxisType>;


template<typename InputType, typename ZAxisType>
using HistMapData = std::tuple<cafemol::ValueWithIndex<InputType>, ZAxisType>;


template<typename InputType, typename OutputType>
class DensityFunction {

public:
	DensityFunction() = default;
	~DensityFunction() = default;

	// -----------------------------------------------------------------------------------
	void set_DataRange(const InputType& data_left_lim, const InputType& data_right_lim){
		if (data_left_lim >= data_right_lim) eout("Invalid range set-up, left >= right");
		data_range[0] = data_left_lim;
		data_range[1] = data_right_lim;
	}

	void set_DataRange(const InputType& data_x_left_lim, const InputType& data_x_right_lim, const InputType& data_y_left_lim, const InputType& data_y_right_lim)  {

		if (data_x_left_lim >= data_x_right_lim) eout("Invalid range set-up, left >= right");
		else if (data_y_left_lim >= data_y_right_lim) eout("Invalid range set-up, left >= right");

		data_xy_range[0] = data_x_left_lim;
		data_xy_range[1] = data_x_right_lim;
		data_xy_range[2] = data_y_left_lim;
		data_xy_range[3] = data_y_right_lim;
	}

	// -----------------------------------------------------------------------------------
	std::vector<HistogramData<InputType, OutputType>> make_Density(std::vector<InputType> vec, const InputType& bin_width) {

		std::sort(vec.begin(), vec.end());
	
		std::vector<HistogramData<InputType, OutputType>> result;
		
		InputType data_left_lim;
		if (data_range.empty()) {
			typename std::vector<InputType>::const_iterator vec_min_itr = std::min_element(vec.begin(), vec.end());
			data_left_lim = *vec_min_itr;
		}
		else {
			data_left_lim = data_range[0];
		}
		
		int count = 0;
		InputType current_point = data_left_lim;
	//	InputType vector_size = vec.size();
	
		for (const InputType& vec_value : vec) {
			if ((current_point <= vec_value) && ((current_point + bin_width) > vec_value)) ++count;
	
			else if (current_point + bin_width <= vec_value) {
				OutputType density = static_cast<OutputType>(count);
				result.push_back(std::make_tuple(current_point, density));
				current_point += bin_width;
				count = 1;
			}
		}
		OutputType density = static_cast<InputType>(count);
		result.push_back(std::make_tuple(current_point, density));
	
		return result;
	}

	// -----------------------------------------------------------------------------------
	std::vector<HistMapData<InputType, OutputType>> make_Density(const std::array<std::vector<InputType>, 2>& vec, const std::array<InputType, 2>& bin_widths) {

		const std::vector<cafemol::ValueWithIndex<InputType>> xy_data = convert_and_sort_XY(vec);
		const std::vector<HistMapData<InputType, OutputType>> result;

		InputType ini_subsection_val;
		if (data_xy_range.empty()) {
			ini_subsection_val = xy_data[0].Index();
		}
		else {
			ini_subsection_val = data_xy_range[0];
		}

		for (const cafemol::ValueWithIndex<InputType>& xy : xy_data) {
			std::vector<InputType> y_data;
			if ((ini_subsection_val <= xy.Index()) && ((ini_subsection_val + bin_widths[0]) > xy.Index())) y_data.push_back(xy.Value());
			else if (ini_subsection_val + bin_widths[0] <= xy.Index()) {
				std::vector<HistogramData<InputType, OutputType>> histogram_on_yaxis 
					= make_Density(y_data, bin_widths[1]);
				for (std::size_t idx = 0; idx < histogram_on_yaxis.size(); ++idx) {
					result.push_back(std::make_tuple(ini_subsection_val,
													 std::get<0>(histogram_on_yaxis[idx]),
													 std::get<1>(histogram_on_yaxis[idx])));
				}
			}
		}
	}


private:

	// -----------------------------------------------------------------------------------
	std::array<InputType, 2> data_range;
	std::array<InputType, 4> data_xy_range;

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();

	// -----------------------------------------------------------------------------------

	std::vector<cafemol::ValueWithIndex<InputType>> convert_and_sort_XY(const std::array<std::vector<InputType>, 2>& vec) {

		if (vec[0].size() != vec[1].size()) eout("the degree of freedom of x-y is not consistent");
		const std::size_t& vector_size = vec[0].size();

		std::vector<cafemol::ValueWithIndex<InputType>> vector_with_indices;
		for (std::size_t idx = 0; idx < vector_size; ++idx) {
			vector_with_indices.push_back(cafemol::ValueWithIndex(vec[0][idx], vec[1][idx]));
		}

		std::sort(vector_with_indices.begin(), vector_with_indices.end(), cafemol::compare_VWI_Index);

		typename std::vector<cafemol::ValueWithIndex<InputType>>::iterator sort_begin_itr =
			vector_with_indices.begin();
		typename std::vector<cafemol::ValueWithIndex<InputType>>::iterator sort_end_itr =
			vector_with_indices.begin();

		for (typename std::vector<cafemol::ValueWithIndex<InputType>>::iterator itr = vector_with_indices.begin(); itr != vector_with_indices.end(); ++itr) {
			if (*itr.Index() == *sort_begin_itr.Index()) continue;

			std::sort(sort_begin_itr, itr, compare_VWI_Value);
			sort_begin_itr = itr;
		}
		std::sort(sort_begin_itr, vector_with_indices.end(), compare_VWI_Value);

		return vector_with_indices;
	}

};

}

#endif // DENSITY_FUNCTION_HPP
