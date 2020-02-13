#include<Eigen/Core>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<OtherFormat/MatrixFileReader.hpp>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<memory>


int main(int argc, char* argv[]) {

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != 6) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_prefix = argv[1];
	const std::size_t& filenum = std::stoi(argv[2]);
	const std::string& s_filenum = argv[2];
	const std::string& input_suffix = argv[3];
	const std::string& output_name = argv[4];
	const double& data_size_in_one = std::stof(argv[5]);
	const std::size_t& filenum_digit = s_filenum.size();

	std::ifstream skip_file("SkipNum", std::ios::in);
	std::vector<std::size_t> skip_number;

	if (skip_file.is_open()) {
		std::string buffer;
		while (std::getline(skip_file, buffer)) {
			skip_number.push_back(std::stoi(buffer));
		}
	}
	skip_file.close();

	Eigen::MatrixXd result;

	std::unique_ptr<cafemol::MatrixFileReader> matrix_reader = std::make_unique<cafemol::MatrixFileReader>();

	bool is_allocated = false;
	std::size_t read_num = 0;
	sout("Calculation starts.");

	for (std::size_t file_index = 1; file_index <= filenum; ++file_index) {
		bool is_skipped = false;

		for (const std::size_t& skip_idx : skip_number) {
			is_skipped = (is_skipped || skip_idx == file_index);
		}
		if (is_skipped) {
			sout("File Index: " + std::to_string(file_index) + " skipped.");
			continue;
		}

		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(filenum_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;
		
		std::string input_name(buffer_filename.str());
		if (!is_allocated) {
			const Eigen::MatrixXd& matrix = matrix_reader->read_Matrix(input_name);
			result.resize(matrix.rows(), matrix.cols());
			result = matrix;
			is_allocated = true;
		}
		else {
			const Eigen::MatrixXd& matrix = matrix_reader->read_Matrix(input_name);
			result += matrix;
		}
		++read_num;
	}

	const double& f_read_num = static_cast<double>(read_num);
	const double& full_data_size = f_read_num * data_size_in_one;

	result /= full_data_size;

	sout("Output to " + output_name);

	std::ofstream ofs(output_name, std::ios::out | std::ios::binary);
	std::int32_t block_size = sizeof(int) * 2;
	int row_size = result.rows();
	int col_size = result.cols();
	ofs.write(reinterpret_cast<char*>(&block_size), sizeof(std::int32_t));
	ofs.write(reinterpret_cast<char*>(&row_size), sizeof(int));
	ofs.write(reinterpret_cast<char*>(&col_size), sizeof(int));
	ofs.write(reinterpret_cast<char*>(&block_size), sizeof(std::int32_t));

	std::int32_t mat_data_size = sizeof(double) * result.size();

	ofs.write(reinterpret_cast<char*>(&mat_data_size), sizeof(std::int32_t));

	for (int i_mat_datum = 0; i_mat_datum < result.size(); ++i_mat_datum) {
		ofs.write(reinterpret_cast<char*>(&result(i_mat_datum)), sizeof(double));
	}

	ofs.write(reinterpret_cast<char*>(&mat_data_size), sizeof(std::int32_t));

	ofs.close();


	return 0;
}
