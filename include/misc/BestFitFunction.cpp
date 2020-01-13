#include"BestFitFunction.hpp"


// public

Eigen::Matrix<float, Eigen::Dynamic, 3> cafemol::library::Best_Fit_Performer::perform_EigenBestFit(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step) {

	// convert STL matrix to Eigen Matrix
	Eigen::Matrix<float, Eigen::Dynamic, 3>&& ini_xyz_mat = cafemol::library::convert_2EigenMat(xyz_ref);
	Eigen::Matrix<float, Eigen::Dynamic, 3>&& xyz_mat = cafemol::library::convert_2EigenMat(xyz_at_each_step);
	
	Eigen::Vector3f ini_xyz_centroid = ini_xyz_mat.colwise().sum() / ini_xyz_mat.rows();
//	ini_xyz_mat.rowwise() -= ini_xyz_centroid;
	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < ini_xyz_mat.rows(); ++idx) {
			ini_xyz_mat(idx, idim) -= ini_xyz_centroid(idim);
		}
	}

	Eigen::Vector3f xyz_centroid = xyz_mat.colwise().sum() / xyz_mat.rows();
//	xyz_mat.rowwise() -= xyz_centroid;
	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < xyz_mat.rows(); ++idx) {
			xyz_mat(idx, idim) -= xyz_centroid(idim);
		}
	}

	Eigen::Matrix3f covariance_mat = xyz_mat.transpose() * ini_xyz_mat;

	Eigen::JacobiSVD<Eigen::Matrix3f> SVD(covariance_mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3f matrix_U = SVD.matrixU();
	Eigen::Matrix3f matrix_V = SVD.matrixV();
	float VU_determinant = (matrix_V * matrix_U.transpose()).determinant();
	Eigen::Matrix3f matrix_S;
	matrix_S << 1.0, 0.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 0.0, VU_determinant;
  
	Eigen::Matrix3f rotation_mat = (matrix_V * matrix_S * matrix_U.transpose()).transpose();

//	Eigen::Matrix3f rotation_mat = (matrix_V * matrix_U.transpose()).transpose();
//
//	if (rotation_mat.determinant() < 0) {
//		matrix_V.row(2) *= -1;
//		rotation_mat = (matrix_V * matrix_U.transpose()).transpose();
//	}

	Eigen::Matrix<float, Eigen::Dynamic, 3>&& xyz_best_fit = xyz_mat * rotation_mat;
//	std::array<std::vector<float>, 3>&& result = cafemol::library::convert_2STLMat(xyz_best_fit);
//	return result;
//	move to the center of a reference structure 
	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < xyz_best_fit.rows(); ++idx) {
			xyz_best_fit(idx, idim) += xyz_centroid(idim);
		}
	}

	return xyz_best_fit;

}


// protected

std::array<std::vector<float>, 3> cafemol::library::Best_Fit_Performer::perform_BestFit(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step) {

	std::array<std::vector<float>, 3> xyz_target;
	for (std::size_t idx = 0; idx < 3; ++idx) {
		xyz_target[idx].resize(xyz_at_each_step[idx].size());
	}
//	std::cout << "Target trajectory size ->"
//			  << xyz_target[0].size() << " "
//			  << xyz_target[1].size() << " "
//			  << xyz_target[2].size() << std::endl;

	// calculation of the center of molecules
	std::array<float, 3> center_of_ref = calc_CenterOfMass(xyz_ref);
	std::array<float, 3> center_at_each_step = calc_CenterOfMass(xyz_at_each_step);
//	std::array<float, 3> translation_vector = center_of_ref - center_at_each_step;

	// calculation of covariance matrix
	Eigen::Matrix<float, 3, 3> covariance_matrix = calc_CrossCovarianceMatrix(xyz_ref, xyz_at_each_step);
	// calculation of rotation matrix
	Eigen::Matrix<float, 3, 3> rotation_matrix = calc_BestFitRotationMatrix(covariance_matrix);

	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < xyz_at_each_step[idim].size(); ++idx) {
			xyz_target[idim][idx] = 0.0;
			for (std::size_t idim_k = 0; idim_k < 3; ++idim_k) {
		//		xyz_target[idim][idx] += rotation_matrix(idim, idim_k) * 
		//								 (xyz_at_each_step[idim_k][idx] + translation_vector[idim_k]);
				xyz_target[idim][idx] += rotation_matrix(idim, idim_k) * 
										 (xyz_at_each_step[idim_k][idx] - center_at_each_step[idim_k]);
			}
			// move to the center of a reference structure
			xyz_target[idim][idx] += center_of_ref[idim];
		}
	}
	
	return xyz_target;
}


Eigen::VectorXd cafemol::library::Best_Fit_Performer::perform_BestFit_ravel_Structure(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step) {

	Eigen::VectorXd xyz_target(3 * xyz_at_each_step[0].size());
//	std::cout << "Target trajectory size ->"
//			  << xyz_target[0].size() << " "
//			  << xyz_target[1].size() << " "
//			  << xyz_target[2].size() << std::endl;

	// calculation of the center of molecules
	std::array<float, 3> center_of_ref = calc_CenterOfMass(xyz_ref);
	std::array<float, 3> center_at_each_step = calc_CenterOfMass(xyz_at_each_step);
//	std::array<float, 3> translation_vector = center_of_ref - center_at_each_step;

	// calculation of covariance matrix
	Eigen::Matrix<float, 3, 3> covariance_matrix = calc_CrossCovarianceMatrix(xyz_ref, xyz_at_each_step);
	// calculation of rotation matrix
	Eigen::Matrix<float, 3, 3> rotation_matrix = calc_BestFitRotationMatrix(covariance_matrix);

	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < xyz_at_each_step[idim].size(); ++idx) {
			xyz_target(idim + 3 * idx) = 0;
			for (std::size_t idim_k = 0; idim_k < 3; ++idim_k) {
		//		xyz_target[idim][idx] += rotation_matrix(idim, idim_k) * 
		//								 (xyz_at_each_step[idim_k][idx] + translation_vector[idim_k]);
				xyz_target(idim + 3 * idx) += rotation_matrix(idim, idim_k) * 
										 (xyz_at_each_step[idim_k][idx] - center_at_each_step[idim_k]);
			}
			// move to the origin
			xyz_target(idim + 3 * idx) -= center_at_each_step[idim];
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


Eigen::Matrix<float, 3, 3> cafemol::library::Best_Fit_Performer::calc_CrossCovarianceMatrix(const std::array<std::vector<float>, 3>& xyz_ref, const std::array<std::vector<float>, 3>& xyz_at_each_step) {

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
				result(idim, jdim) += (xyz_ref[jdim][idx] - center_of_ref[jdim]) *
									  (xyz_at_each_step[idim][idx] - center_at_each_step[idim]);
			}
		}
	}

	return result;

}


Eigen::Matrix<float, 3, 3> cafemol::library::Best_Fit_Performer::calc_BestFitRotationMatrix(const Eigen::Matrix<float, 3, 3>& covariance_matrix) {
	// perform SVD
	Eigen::JacobiSVD<Eigen::Matrix<float, 3, 3>> SVD(covariance_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix<float, 3, 3> matrix_U = SVD.matrixU();
	Eigen::Matrix<float, 3, 3> matrix_V = SVD.matrixV();
//	Eigen::Matrix<float, 3, 3> rotation_matrix = matrix_V * matrix_U.transpose();

	float det_VU = (matrix_V * matrix_U.transpose()).determinant();
	Eigen::Matrix<float, 3, 3> matrix_S;
	matrix_S << 1.0, 0.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 0.0, det_VU;
	Eigen::Matrix<float, 3, 3> rotation_matrix = matrix_V * matrix_S * matrix_U.transpose();


	// adjust the mirror image
//	if (rotation_matrix.determinant() < 0) {
//		matrix_V.row(2) *= -1;
//	//	rotation_matrix = (matrix_V * matrix_U.transpose()).transpose();
//		rotation_matrix = matrix_V * matrix_U.transpose();
//	}

	return rotation_matrix;
}
