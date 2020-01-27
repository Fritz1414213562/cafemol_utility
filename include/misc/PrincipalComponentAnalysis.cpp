#include"PrincipalComponentAnalysis.hpp"


// private


std::array<std::vector<float>, 3> cafemol::analysis::PCA_Performer::calc_AverageStructure(const cafemol::Trajectory& traj) {
//std::array<std::vector<float>, 3> cafemol::library::PCA_Performer::calc_AverageStructure(const cafemol::Trajectory& traj) {

	std::array<std::vector<float>, 3> average_structure;
	for (std::size_t idim = 0; idim < 3; ++idim) {
		average_structure[idim].resize(traj[0][idim].size(), 0.0);
	}

	std::size_t n_frame = 0;
	for (const std::array<std::vector<float>, 3>& snap_shot : traj) {
		for (std::size_t idim = 0; idim < 3; ++idim) {
			for (std::size_t idx = 0; idx < snap_shot[idim].size(); ++idx) {
				average_structure[idim][idx] += snap_shot[idim][idx];
			}
		}
		++n_frame;
	}

	for (std::size_t idim = 0; idim < 3; ++idim) {
		for (std::size_t idx = 0; idx < average_structure[idim].size(); ++idx) {
			average_structure[idim][idx] /= n_frame;
		}
	}

	return average_structure;
}


Eigen::VectorXd cafemol::analysis::PCA_Performer::convert_and_ravel_XYZStructure(const std::array<std::vector<float>, 3>& xyz) {
//Eigen::VectorXd cafemol::library::PCA_Performer::convert_and_ravel_XYZStructure(const std::array<std::vector<float>, 3>& xyz) {

	if ((xyz[0].size() != xyz[1].size()) || (xyz[0].size() != xyz[2].size())) eout("The degree of freedom of xyz must be the same.");

	const std::size_t& atom_number = xyz[0].size();
	Eigen::VectorXd result(3 * atom_number);
	const std::array<float, 3>& center_of_mass = calc_CenterOfMass(xyz);
	for (std::size_t i_atom = 0; i_atom < atom_number; ++i_atom) {
		for (std::size_t i_dim = 0; i_dim < 3; ++i_dim) {
			const std::size_t& idx = 3 * i_atom + i_dim;
			result(idx) = static_cast<double>(xyz[i_dim][i_atom] - center_of_mass[i_dim]);
		}
	}

//	if (is_moved2origin) {
//		const std::array<float, 3>& center_of_mass = calc_CenterOfMass(xyz);
//		for (std::size_t i_atom = 0; i_atom < atom_number; ++i_atom) {
//			for (std::size_t i_dim = 0; i_dim < 3; ++i_dim) {
//				const std::size_t& idx = 3 * i_atom + i_dim;
//				result(idx) = static_cast<double>(xyz[i_dim][i_atom] - center_of_mass[i_dim]);
//			}
//		}
//	}
//	else {
//		for (std::size_t i_atom = 0; i_atom < atom_number; ++i_atom) {
//			for (std::size_t i_dim = 0; i_dim < 3; ++i_dim) {
//				const std::size_t& idx = 3 * i_atom + i_dim;
//				result(idx) = static_cast<double>(xyz[i_dim][i_atom]);
//			}
//		}
//	}

	return result;
}


Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cafemol::analysis::PCA_Performer::calc_CovarianceMatrix(const std::array<std::vector<float>, 3>& xyz) {
//Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cafemol::library::PCA_Performer::calc_CovarianceMatrix(const std::array<std::vector<float>, 3>& xyz) {


	const std::size_t atom_number = xyz[0].size();
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result(3 * atom_number, 3 * atom_number);

	if (is_moved2origin) {
		const std::array<float, 3>&& center_of_mass = calc_CenterOfMass(xyz);
		for (std::size_t i_atom = 0; i_atom < atom_number; ++i_atom) {
			for (std::size_t i_dim = 0; i_dim < 3; ++i_dim) {
				const std::size_t& i_idx = 3 * i_atom + i_dim;

				for (std::size_t j_atom = 0; j_atom < atom_number; ++j_atom) {
					for (std::size_t j_dim = 0; j_dim < 3; ++j_dim) {
						const std::size_t& j_idx = 3 * j_atom + j_dim;

						result(i_idx, j_idx) = static_cast<double>(
											   xyz[i_dim][i_atom] - center_of_mass[i_dim]) *
											   static_cast<double>(
											   xyz[j_dim][j_atom] - center_of_mass[j_dim]);
					}
				}
			}
		}
	}
	else {
		for (std::size_t i_atom = 0; i_atom < atom_number; ++i_atom) {
			for (std::size_t i_dim = 0; i_dim < 3; ++i_dim) {
				const std::size_t& i_idx = 3 * i_atom + i_dim;

				for (std::size_t j_atom = 0; j_atom < atom_number; ++j_atom) {
					for (std::size_t j_dim = 0; j_dim < 3; ++j_dim) {
						const std::size_t& j_idx = 3 * j_atom + j_dim;

						result(i_idx, j_idx) = static_cast<double>(
											   xyz[i_dim][i_atom]) *
											   static_cast<double>(
											   xyz[j_dim][j_atom]);
					}
				}
			}
		}
	}

	return result;
}


// You must ensure that each snapshot of the trajectory file have the same size among xyz directions.
cafemol::PrincipalComponents cafemol::analysis::PCA_Performer::perform_PrincipalComponentAnalysis(const cafemol::Trajectory& traj, const std::size_t& max_component_num) {
//cafemol::PrincipalComponents cafemol::library::PCA_Performer::perform_PrincipalComponentAnalysis(const cafemol::Trajectory& traj, const std::size_t& max_component_num) {


	const std::size_t& frame_number = traj.size();
	if (frame_number <= 1) eout("The frame size of a trajectory must be more than 1.");
	sout("The frame size except the frame on the reference structure is " + std::to_string(frame_number) + ".",
		"");

	const std::size_t& atom_number = traj[0][0].size();
	if (max_component_num > 3 * atom_number) eout("The atom size is smaller than component size.");

	sout("The atom size is " + std::to_string(atom_number) + ".", "");

	// prepare results
	cafemol::PrincipalComponents result(max_component_num);

	// save the reference structure
	const std::array<std::vector<float>, 3>& xyz_ref = traj[0];
	// save the best-fit structures
	cafemol::EigenTrajectory best_fit_traj;

	// calculate the covariance matrices of all snapshots and take their average.
	sout("Calculating the average Covariance matrix");
	Eigen::MatrixXd average_covariance_matrix(3 * atom_number, 3 * atom_number);
	for (std::size_t iframe = 1; iframe < frame_number; ++iframe) {
	//	if (iframe % calculate_step != 0) continue;
		const std::array<std::vector<float>, 3>& xyz = traj[iframe];
		const std::array<std::vector<float>, 3>& xyz_fit_on_ref = perform_BestFit(xyz_ref, xyz);
		const Eigen::VectorXd& xyz_raveled = convert_and_ravel_XYZStructure(xyz_fit_on_ref);
		best_fit_traj.push_back(xyz_raveled);
	//	average_covariance_matrix += calc_CovarianceMatrix(xyz_fit_on_ref);
		average_covariance_matrix += xyz_raveled * xyz_raveled.transpose();
	}
//	average_covariance_matrix /= ((frame_number - 1) / calculate_step);
	average_covariance_matrix /= (frame_number - 1);
	sout(".... Done");

	sout("solving the eigen values and vectors");
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(average_covariance_matrix);
	// the sum of contribution rates
	const double& contribution_rates_sum = eigen_solver.eigenvalues().sum();

	for (std::size_t i_component = 0; i_component < max_component_num; ++i_component) {
		const std::size_t& eigen_value_idx = 3 * atom_number - (1 + i_component);
		result.contribution_rates(i_component) = eigen_solver.eigenvalues()(eigen_value_idx) / 
												 contribution_rates_sum;

		result.principal_component_axis[i_component] = eigen_solver.eigenvectors().col(eigen_value_idx);
	}
	sout(".... Done", "");

  // project the best-fit trajectory on PCs;

	std::size_t skip_iframe = 0;
	for (std::size_t iframe = 0; iframe < frame_number - 1; ++iframe) {
  	// calculation methods for stl container
		if (iframe % calculate_step != 0) continue;
  //	const std::array<std::vector<float>, 3>& xyz_fit_on_ref = best_fit_traj[skip_iframe];
  //	const std::array<std::vector<float>, 3>& xyz_fit_on_ref = best_fit_traj[iframe];
  //	result.projections.push_back(project_StructureOnPCs(xyz_fit_on_ref, result.principal_component_axis));
		++skip_iframe;
  	
  	// calculation methods for eigen
		const Eigen::VectorXd& xyz_fit_on_ref = best_fit_traj[iframe];
		std::vector<double> projected_position;
		for (std::size_t i_component = 0; i_component < max_component_num; ++i_component) {
			const Eigen::VectorXd& pc_vector = result.principal_component_axis[i_component];
			const double& projected_coordinate = pc_vector.transpose() * xyz_fit_on_ref;
			projected_position.push_back(projected_coordinate);
		}
		result.projections.push_back(projected_position);
	}

	return result;

}

cafemol::PrincipalComponents cafemol::analysis::PCA_Performer::perform_PrincipalComponentAnalysis(const cafemol::EigenTrajectory& traj, const std::size_t& max_component_num) {
//cafemol::PrincipalComponents cafemol::library::PCA_Performer::perform_PrincipalComponentAnalysis(const cafemol::EigenTrajectory& traj, const std::size_t& max_component_num) {
// for EigenTrajectory (You must ensure that the trajectory has been fit on the reference)

	const std::size_t& atom_number = traj[0].size();
	if (max_component_num > atom_number) eout("The atom size is smaller than component size.");

	sout("The degree of freedom of this system is " + std::to_string(atom_number) + ".", "");

	// prepare result
	cafemol::PrincipalComponents result(max_component_num);

	const std::size_t& frame_number = traj.size();

	// calculate the covariance matrices of all snapshots and take their average.
	sout("Calculating the average covariance matrix");
	Eigen::MatrixXd average_covariance_matrix(atom_number, atom_number);

	for (std::size_t iframe = 0; iframe < frame_number; ++iframe) {
		if (iframe % calculate_step != 0) continue;
//		if (iframe % 1000 == 0) sout(iframe);
		average_covariance_matrix += traj[iframe] * traj[iframe].transpose();
	}
	average_covariance_matrix /= (frame_number / calculate_step);
	sout(".... Done");

	sout("solving the eigen values and vectors");
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(average_covariance_matrix);
	// the sum of contribution rates
	const double& contribution_rates_sum = eigen_solver.eigenvalues().sum();

	for (std::size_t i_component = 0; i_component < max_component_num; ++i_component) {
		const std::size_t& eigen_value_idx = atom_number - (1 + i_component);
		result.contribution_rates(i_component) = eigen_solver.eigenvalues()(eigen_value_idx) / 
												 contribution_rates_sum;
		result.principal_component_axis[i_component] = eigen_solver.eigenvectors().col(eigen_value_idx);
	}
	sout(".... Done");

	// project the best-fit trajectory on PCs

	for (std::size_t iframe = 0; iframe < frame_number; ++iframe) {
		if (iframe % calculate_step != 0) continue;
		const Eigen::VectorXd& xyz_fit_on_ref = traj[iframe];
		std::vector<double> projected_position;
		for (std::size_t i_component = 0; i_component < max_component_num; ++i_component) {
			const Eigen::VectorXd& pc_vector = result.principal_component_axis[i_component];
			const double& projected_coordinate = pc_vector.transpose() * xyz_fit_on_ref;
			projected_position.push_back(projected_coordinate);
		}
		result.projections.push_back(projected_position);
	}

	return result;

}


std::vector<double> cafemol::analysis::PCA_Performer::project_StructureOnPCs(const std::array<std::vector<float>, 3>& xyz, const std::vector<Eigen::VectorXd>& PC_axis) {
//std::vector<double> cafemol::library::PCA_Performer::project_StructureOnPCs(const std::array<std::vector<float>, 3>& xyz, const std::vector<Eigen::VectorXd>& PC_axis) {

	std::vector<double> result(PC_axis.size(), 0);

	if (is_moved2origin) {
		const std::size_t& atom_size = xyz[0].size();
		const std::array<float, 3>&& center_of_mass = calc_CenterOfMass(xyz);
		for (std::size_t i_atom = 0; i_atom < atom_size; ++i_atom) {
			for (std::size_t idim = 0; idim < 3; ++idim) {
				const std::size_t& idx = 3 * i_atom + idim;
				for (std::size_t i_pc = 0; i_pc < PC_axis.size(); ++i_pc) {
					result[i_pc] += static_cast<double>(xyz[idim][i_atom] - center_of_mass[idim]) * 
														PC_axis[i_pc](idx);
				}
			}
		}
	}
	else {
		const std::size_t& atom_size = xyz[0].size();
		for (std::size_t i_atom = 0; i_atom < atom_size; ++i_atom) {
			for (std::size_t idim = 0; idim < 3; ++idim) {
				const std::size_t& idx = 3 * i_atom + idim;
				for (std::size_t i_pc = 0; i_pc < PC_axis.size(); ++i_pc) {
					result[i_pc] += static_cast<double>(xyz[idim][i_atom]) * PC_axis[i_pc](idx);
				}
			}
		}
	}

	return result;
}



void cafemol::analysis::PCA_Performer::run(const cafemol::Trajectory& traj, const std::size_t& max_component_num) {
//void cafemol::library::PCA_Performer::run(const cafemol::Trajectory& traj, const std::size_t& max_component_num) {

	// open or make output file
	if (output_name.empty()) eout("You must initialize output file.");

	std::ofstream ofs(output_name, std::ios::out);

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Principal Components Analysis starts.", "");
	sout("Best fit structures to the reference");

	// ravel and best fit the structure on the reference
	cafemol::EigenTrajectory best_fit_traj;
	const std::size_t& frame_number = traj.size();

	if (frame_number <= 1) eout("The frame size of a trajectory must be more than 1.");
	sout("The frame size except the frame on the reference structure is " + std::to_string(frame_number - 1) + ".",
		"");
	const std::array<std::vector<float>, 3>& xyz_ref = traj[0];

	sout("convert and fit the structure on the reference");
	for (std::size_t iframe = 1; iframe < frame_number; ++iframe) {
		const std::array<std::vector<float>, 3>& xyz = traj[iframe];
		best_fit_traj.push_back(perform_BestFit_ravel_Structure(xyz_ref, xyz));
	}

	cafemol::PrincipalComponents pc_result = perform_PrincipalComponentAnalysis(best_fit_traj, 
																				max_component_num);

	sout.output_HyphenBlock("output the result to " + output_name, BLOCK_SIZE);

	ofs << "<<<< n_principal_components" << std::endl << pc_result.size() << std::endl;
	ofs << ">>>>" << std::endl << std::endl;

	ofs << "<<<< contribution_rates" << std::endl;
	for (std::size_t i_axis = 0; i_axis < pc_result.size(); ++i_axis) {
		ofs << "PC" << i_axis + 1 << ": " << pc_result.contribution_rates(i_axis) << std::endl;
	}
	ofs << ">>>>" << std::endl << std::endl;

	ofs << "<<<< principal_component_vectors" << std::endl;
	ofs << "degree of freedom, PC1 vector, PC2 vector" << std::endl;
	const std::size_t& pc_vec_size = pc_result.principal_component_axis[0].size();
	for (std::size_t idx = 0; idx < pc_vec_size; ++idx) {
		ofs << idx + 1;
		for (std::size_t i_axis = 0; i_axis < pc_result.size(); ++i_axis) {
			ofs << " " << pc_result.principal_component_axis[i_axis](idx);
		}
		ofs << std::endl;
	}
	ofs << ">>>>" << std::endl << std::endl;

	ofs << "<<<< projection_on_PCs" << std::endl;
	ofs << "MD frame, PC1, PC2" << std::endl;
	for (std::size_t iframe = 0; iframe < pc_result.projections.size(); ++iframe) {
//		ofs << iframe + 1;
		ofs << (iframe + 1) * calculate_step;
		const std::vector<double>& projection = pc_result.projections[iframe];
		for (std::size_t idx = 0; idx < projection.size(); ++idx) {
			ofs << " " << projection[idx];
		}
		ofs << std::endl;
	}
	ofs << ">>>>" << std::endl;

	sout(".... Done");
	ofs.close();

}

//void cafemol::library::PCA_Performer::run(const cafemol::EigenTrajectory& traj, const std::size_t& max_component_num) {
//
//	// open or make output file
//	if (output_name.empty()) eout("You must initialize output file.");
//
//	std::ofstream ofs(output_name, std::ios::out);
//
//	sout.output_HyphenBlock("", BLOCK_SIZE);
//	sout("Principal Components Analysis starts.", "");
//
//	cafemol::PrincipalComponents pc_result = perform_PrincipalComponentAnalysis(traj, 
//																				max_component_num);
//
//	sout.output_HyphenBlock("output the result to " + output_name, BLOCK_SIZE);
//
//	ofs << "<<<< the number of principal components" << std::endl << pc_result.size() << std::endl;
//	ofs << ">>>>" << std::endl << std::endl;
//
//	ofs << "<<<< the contribution rates of principal components" << std::endl;
//	for (std::size_t i_axis = 0; i_axis < pc_result.size(); ++i_axis) {
//		ofs << "PC" << i_axis + 1 << ": " << pc_result.contribution_rates(i_axis);
//	}
//	ofs << ">>>>" << std::endl << std::endl;
//
//	ofs << "<<<< the principal components vectors" << std::endl;
//	ofs << "degree of freedom, PC1 vector, PC2 vector" << std::endl;
//	const std::size_t& pc_vec_size = pc_result.principal_component_axis[0].size();
//	for (std::size_t idx = 0; idx < pc_vec_size; ++idx) {
//		ofs << idx + 1;
//		for (std::size_t i_axis = 0; i_axis < pc_result.size(); ++i_axis) {
//			ofs << " " << pc_result.principal_component_axis[i_axis](idx);
//		}
//		ofs << std::endl;
//	}
//	ofs << ">>>>" << std::endl << std::endl;
//
//	ofs << "<<<< the structure projection on PC plane at each MD frame" << std::endl;
//	ofs << "MD frame, PC1, PC2" << std::endl;
//	for (std::size_t iframe = 0; iframe < pc_result.projections.size(); ++iframe) {
////		ofs << (iframe + 1) * skip_iframe;
//		ofs << iframe + 1;
//		const std::vector<double>& projection = pc_result.projections[iframe];
//		for (std::size_t idx = 0; idx < projection.size(); ++idx) {
//			ofs << " " << projection[idx];
//		}
//		ofs << std::endl;
//	}
//	ofs << ">>>>" << std::endl;
//
//	sout(".... Done");
//	ofs.close();
//
//}
