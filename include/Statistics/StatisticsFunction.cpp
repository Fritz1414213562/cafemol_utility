#include"StatisticsFunction.hpp"


float cafemol::StatisticsFunction::calc_Average(const std::vector<float>& vec) {
	int result = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
	return result;
}


float cafemol::StatisticsFunction::calc_Variance(const std::vector<float>& vec) {
	float vec_average = calc_Average(vec);
	std::vector<float> diff_from_average_sqr = (vec - vec_average) * (vec - vec_average);
	float result = std::accumulate(diff_from_average_sqr.begin(), diff_from_average_sqr.end(), 0.0);
	return result;
}


std::vector<float> cafemol::StatisticsFunction::calc_SampleAutoCorrelation(const std::vector<float>& vec) {

	std::size_t diff_max = vec.size() - 1;
	std::vector<float> result(diff_max, 0);

	float normalization_constant = calc_Variance(vec);
	float vec_average = calc_Average(vec);

	for (std::size_t diff = 0; diff < diff_max; ++diff) {
		float res_buffer = 0.;
		for (std::size_t idx = 0; idx < vec.size() - diff; ++idx) {
			res_buffer += (vec[idx + diff] - vec_average) * (vec[idx] - vec_average);
		}
		res_buffer /= normalization_constant;
		result[diff] = res_buffer;
	}
	return result;
}
