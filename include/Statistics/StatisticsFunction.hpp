#ifndef STATISTICS_FUNCTION_HPP
#define STATISTICS_FUNCTION_HPP

#include<iostream>
#include<fstream>
#include<vector>
#include<numeric>
#include<algorithm>

namespace cafemol {

class StatisticsFunction {

public:
	StatisticsFunction() = default;
	~StatisticsFunction() = default;
	std::vector<float> calc_SampleAutoCorrelation(const std::vector<float>& vec);

protected:
	float calc_Average(const std::vector<float>& vec);
	float calc_Variance(const std::vector<float>& vec);

};


// inline function

template<typename realT>
inline std::vector<realT> operator-(const std::vector<realT>& lhs, const realT& rhs) {
	std::vector<realT> result(lhs.size(), 0);
	for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
		result[idx] = lhs[idx] - rhs;
	}
	return result;
}


template<typename realT>
inline std::vector<realT> operator*(const std::vector<realT>& lhs, const std::vector<realT>& rhs) {

	try {
		bool is_consistent = (lhs.size() == rhs.size());
		if (!is_consistent) {
			throw std::exception();
		}
	} catch (std::exception e) {
		std::cerr <<"Error: The vector size is not consistent with another one." << std::endl;
		std::exit(1);
	}
	std::vector<realT> result(lhs.size(), 0);
	for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
		result[idx] = lhs[idx] * rhs[idx];
	}
	return result;
}

}

#endif /* STATISTICS_FUNCTION_HPP */
