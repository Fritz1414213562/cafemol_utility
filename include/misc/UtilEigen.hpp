#ifndef UTIL_EIGEN_HPP
#define UTIL_EIGEN_HPP
#include<Eigen/Core>
#include<iostream>
#include<array>
#include<vector>

namespace cafemol {

namespace library {

template<typename realT, std::size_t column_size>
inline Eigen::Matrix<realT, Eigen::Dynamic, column_size> convert_2EigenMat(const std::array<std::vector<realT>, column_size>& xyz) {

	if ((xyz[0].size() != xyz[1].size()) && (xyz[0].size() != xyz[2].size())) {
		std::cerr << "Error: The arrays shape is not matrix" << std::endl;
		std::exit(1);
	}
	std::size_t row_size = xyz[0].size();
	Eigen::Matrix<realT, Eigen::Dynamic, column_size> result(row_size, column_size);
	for (std::size_t idx = 0; idx < column_size; ++idx) {
		for (std::size_t jdx = 0; jdx < row_size; ++jdx) {
			result(jdx, idx) = xyz[idx][jdx];
		}
	}
	
	return result;
}

template<typename realT, std::size_t column_size>
inline std::array<std::vector<realT>, column_size> convert_2STLMat(const Eigen::Matrix<realT, Eigen::Dynamic, column_size>& xyz) {

	std::array<std::vector<realT>, column_size> result;
	for (std::size_t idx = 0; idx < column_size; ++idx) {
		for (std::size_t jdx = 0; jdx < xyz.rows(); ++jdx) {
			result[idx].push_back(xyz(jdx, idx));
		}
	}

	return result;
}


}
}

#endif /* UTIL_EIGEN_HPP */
