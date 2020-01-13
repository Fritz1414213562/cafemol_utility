#ifndef PRINCIPAL_COMPONENT_ANALYSIS_HPP
#define PRINCIPAL_COMPONENT_ANALYSIS_HPP
#include<BestFitFunction.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<array>
#include<vector>
#include<string>
#include<fstream>


namespace cafemol {

using Trajectories = std::vector<std::array<std::vector<float>, 3>>;
using EigenTrajectory = std::vector<Eigen::VectorXd>;

namespace library {


struct PCA_Container {

public:
	PCA_Container(const std::size_t& axis_size) : axis_number(axis_size) {
		contribution_rates.resize(axis_size);
		principal_component_axis.resize(axis_size);
	//	projections.resize(axis_size)
	}

	Eigen::VectorXd contribution_rates;
	std::vector<Eigen::VectorXd> principal_component_axis;
	std::vector<std::vector<double>> projections;

	const std::size_t& size() {return axis_number;}

private:
	const std::size_t axis_number;

};


class PCA_Performer : public Best_Fit_Performer {

public:
	
	PCA_Performer(const std::string& output_filename) : output_name(output_filename) {}
	~PCA_Performer() = default;

	void set_Trajectory2Origin() {
		is_moved2origin = true;
	}

	void set_CalculationStep(const std::size_t& calc_step) {calculate_step = calc_step;}

	void run(const Trajectories& traj, const std::size_t& max_component_num);
//	void run(const EigenTrajectory& traj, const std::size_t& max_component_num);

private:
	
	// error output
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	// block_size
	const std::size_t BLOCK_SIZE = 80;
	// output name
	const std::string output_name;

	// calculation flags for covariance matrix
	bool is_moved2origin = false;
	std::size_t calculate_step = 1;

	// private methods

	std::array<std::vector<float>, 3> calc_AverageStructure(const Trajectories& traj);

	Eigen::VectorXd convert_and_ravel_XYZStructure(const std::array<std::vector<float>, 3>& xyz);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> calc_CovarianceMatrix(const std::array<std::vector<float>, 3>& xyz);

	PCA_Container perform_PrincipalComponentAnalysis(const Trajectories& traj, const std::size_t& max_component_num);
	PCA_Container perform_PrincipalComponentAnalysis(const EigenTrajectory& traj, const std::size_t& max_component_num);

	std::vector<double> project_StructureOnPCs(const std::array<std::vector<float>, 3>& xyz, const std::vector<Eigen::VectorXd>& PC_axis);

};
}

using PrincipalComponents = library::PCA_Container;

}

#endif /* PRINCIPAL_COMPONENT_ANALYSIS_HPP */
