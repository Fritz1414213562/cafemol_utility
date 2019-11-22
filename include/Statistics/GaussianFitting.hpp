#ifndef GAUSSIAN_FITTING_HPP
#define GAUSSIAN_FITTING_HPP

#include<Eigen/Core>
#include<Eigen/LU>
#include<vector>
#include<array>
#include<iostream>
#include<numeric>
#include<cmath>


namespace cafemol {

class GaussianFitting {

public:
	GaussianFitting() = default;
	~GaussianFitting() = default;

	std::array<float, 3> fit_Curve(const std::vector<float>& data_x, const std::vector<float>& data_y);
	

private:
	Eigen::Matrix3f fitting_matrix;
	Eigen::Vector3f outputs_vector;
	
	void initialize_FittingMatrix(const std::vector<float>& data_x, const std::vector<float>& data_y);
	void initialize_OutputsVector(const std::vector<float>& data_x, const std::vector<float>& data_y);
	Eigen::Vector3f solve_LinearSystem();

};


template<typename realT>
inline std::vector<realT> operator*(const std::vector<realT>& lhs, const std::vector<realT>& rhs) {
	if (lhs.size() != rhs.size()) {
		std::cerr << "Error: The vector size is not consistent with another vector." << std::endl;
		std::exit(1);
	}
	std::vector<realT> res(lhs.size(), 0);
	for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
		res[idx] = lhs[idx] * rhs[idx];
	}
	return res;
}


}
#endif /* GAUSSIAN_FITTING_HPP */
