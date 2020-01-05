#ifndef BEST_FIT_FUNCTION_HPP
#define BEST_FIT_FUNCTION_HPP
#include"ErrorMessage.hpp"
#include<Eigen/Core>
#include<Eigen/SVD>
#include<array>
#include<vector>


namespace cafemol {

namespace library {

class Best_Fit_Performer {

public:

	Best_Fit_Performer() = default;
	~Best_Fit_Performer() = default;

	std::array<std::vector<float>, 3> operator()(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step) {
		return perform_BestFit(xyz_ref, xyz_at_each_step);
	}

protected:

	std::array<std::vector<float>, 3> perform_BestFit(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step);
	
private:

	std::array<float, 3> calc_CenterOfMass(const std::array<std::vector<float>, 3>& xyz);

	Eigen::Matrix<float, 3, 3> calc_CovarianceMatrix(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step);

	Eigen::Matrix<float, 3, 3> calc_BestFitRotationMatrix(const Eigen::Matrix<float, 3, 3>& covariance_matrix);

	// error output
	cafemol::library::Error_Output eout = cafemol::library::Error_Output();

	template<typename realT, std::size_t dim>
	std::array<realT, dim> operator-(const std::array<realT, dim>& lhs, const std::array<realT, dim>& rhs) {
		std::array<realT, dim> result;
		for (std::size_t idim = 0; idim < dim; ++idim) {
			result[idim] = lhs[idim] - rhs[idim];
		}

		return result;
	}

};
}
}

#endif /* BEST_FIT_FUNCTION_HPP */
