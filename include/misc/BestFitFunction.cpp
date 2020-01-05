#include"BestFitFunction.hpp"


// public

std::array<std::vector<float>, 3> perform_BestFit(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step) {

	std::array<std::vector<float>, 3> xyz_target;
	for (std::size_t idx = 0; idx < 3; ++idx) {
		xyz_target[idx].resize(xyz_at_each_step.size());
	}

	// calculation of the center of molecules
	std::array<float, 3> center_of_ref = calc_CenterOfMass(xyz_ref);
	std::array<float, 3> center_at_each_step = calc_CenterOfMass(xyz_at_each_step);
	std::array<float, 3> translation_vector = center_of_ref - center_at_each_step;

	// calculation of covariance matrix
	Eigen::Matrix<float, 3, 3> covariance_matrix = calc_CovarianceMatrix(xyz_ref, xyz_at_each_step);
	// calculation of rotation matrix
	Eigen::Matrix<float, 3, 3> rotation_matrix = calc_BestFitRotationMatrix(covariance_matrix);

	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < xyz_at_each_step.size(); ++idx) {
			xyz_target[idim][idx] = 0.0;
			for (std::size_t idim_k = 0; idim_k < 3; ++idim_k) {
				xyz_target[idim][idx] += rotation_matrix[idim][idim_k] * 
										 (xyz_at_each_step[idim_k][idx] + translation_vector[idim_k];
			}
		}
	}
	
	return xyz_target;
}


// private


std::array<float, 3> cafemol::library::Best_Fit_Performer::calc_CenterOfMass(const std::array<std::vector<float>, 3>& xyz) {

	std::array<float, 3> center_of_mass = {0.0, 0.0, 0.0};

	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < xyz[idim].size(); ++idx) {
			center_of_mass[idim] += xyz[idim][idx];
		}
		center_of_mass[idim] /= xyz[idim].size();
	}

	return center_of_mass;
}


Eigen::Matrix<float, 3, 3> cafemol::library::Best_Fit_Performer::calc_CovarianceMatrix(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step) {

	const std::array<float, 3> center_of_ref = calc_CenterOfMass(xyz_ref);
	const std::array<float, 3> center_at_each_step = calc_CenterOfMass(xyz_at_each_step);

	Eigen::Matrix<float, 3, 3> result;
	result << 0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0;

	// check the vector size
	const std::size_t vector_size = xyz_ref[0].size();
	for (std::size_t idim = 0; idim < 3; ++idim) {
		if ((vector_size == xyz_ref[idim].size()) && (vector_size == xyz_at_each_step[idim].size())) {}
		else eout("The vector size of each snapshot is not consistent with that of reference structure");
	}

	// calculate covariance matrix
	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t jdim = 0; jdim < 3; ++jdim) {
			for (std::size_t idx = 0; idx < vector_size; ++idx) {
				result[idim][jdim] += (xyz_ref[idim][idx] - center_of_ref[idim]) *
									  (xyz_at_each_step[jdim][idx] - center_at_each_step[jdim]);
			}
		}
	}

	return result;

}


Eigen::Matrix<float, 3, 3> calc_BestFitRotationMatrix(const Eigen::Matrix<float, 3, 3>& covariance_matrix) {
	// perform SVD
	Eigen::JacobiSVD<Eigen::Matrix<float, 3, 3>> SVD(covariance_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix<float, 3, 3> matrix_U = SVD.matrixU();
	Eigen::Matrix<float, 3, 3> matrix_V = SVD.matrixV();
	Eigen::Matrix<float, 3, 3> rotation_matrix = (matrix_V * matrix_U.transpose()).transpose();

	// check mirrar image

	if (rotation_matrix.determinant() < 0) {
		matrix_V.row(2) *= -1;
		rotation_matrix = (matrix_V * matrix_U.transpose()).transpose();
	}

	return rotation_matrix;
}
