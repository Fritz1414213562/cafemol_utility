#include<Eigen/Core>
#include<Eigen/Dense>
#include<OtherFormat/MatrixFileReader.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<string>
#include<fstream>
#include<memory>


int main(int argc, char* argv[]) {

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != 4) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& matrix_name = argv[1];
	const std::string& output_binary_name = argv[2];
	const std::string& output_ascii_name = argv[3];

	sout("Calculation starts.");
	sout("Read " + matrix_name);

	std::unique_ptr<cafemol::MatrixFileReader> matrix_reader = std::make_unique<cafemol::MatrixFileReader>();

	const Eigen::MatrixXd& cross_cov = matrix_reader->read_Matrix(matrix_name);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(cross_cov);
	sout(".... Done");

	const double& contribution_rates_sum = eigen_solver.eigenvalues().sum();
	const Eigen::VectorXd& contribution_rates = eigen_solver.eigenvalues() / contribution_rates_sum;
	Eigen::MatrixXd pca_axes = eigen_solver.eigenvectors();

	std::ofstream ascii_file(output_ascii_name, std::ios::out);
	ascii_file << "contribution_rates_sum " << contribution_rates_sum << std::endl;
	ascii_file << "contribution_rates" << std::endl;
	ascii_file << contribution_rates << std::endl;
	ascii_file.close();

	std::ofstream binary_file(output_binary_name, std::ios::out | std::ios::binary);
	std::int32_t block_size = sizeof(int) * 2;
	int row_size = pca_axes.rows();
	int col_size = pca_axes.cols();
	binary_file.write(reinterpret_cast<char*>(&block_size), sizeof(std::int32_t));
	binary_file.write(reinterpret_cast<char*>(&row_size), sizeof(int));
	binary_file.write(reinterpret_cast<char*>(&col_size), sizeof(int));
	binary_file.write(reinterpret_cast<char*>(&block_size), sizeof(std::int32_t));

	std::int32_t mat_data_size = sizeof(double) * pca_axes.size();

	binary_file.write(reinterpret_cast<char*>(&mat_data_size), sizeof(std::int32_t));

	for (int i_mat_datum = 0; i_mat_datum < pca_axes.size(); ++i_mat_datum) {
		binary_file.write(reinterpret_cast<char*>(&pca_axes(i_mat_datum)), sizeof(double));
	}

	binary_file.write(reinterpret_cast<char*>(&mat_data_size), sizeof(std::int32_t));

	binary_file.close();


	return 0;
}
