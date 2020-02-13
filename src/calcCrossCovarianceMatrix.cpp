#include<OtherFormat/ContactStateReader.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<UtilFunc.hpp>
#include<Eigen/Core>
#include<iostream>
#include<fstream>
#include<vector>


int main(int argc, char* argv[]) {

	if (argc != 7) {
		std::cerr << "too much or less arguments" << std::endl;
	}

	const std::string& input_name = argv[1];
	const std::string& average_name = argv[2];
	const std::string& output_name = argv[3];
	const std::size_t& frame_size = std::stoi(argv[4]);
	const std::size_t& data_frame_x = std::stoi(argv[5]);
	const std::size_t& data_frame_y = std::stoi(argv[6]);

	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();


	// ----------------------------------------------------------------------
	// read input

	std::ifstream input_file(input_name, std::ios::in);

	std::string buffer;
	bool is_reading_1stFrame = false;
	Eigen::VectorXd&& data_on_frame = Eigen::VectorXd::Zero(data_frame_x * data_frame_y);
	Eigen::MatrixXd all_data(frame_size, data_frame_x * data_frame_y);

	sout("Reading Contact state data");
	
	std::size_t i_frame = 0;
	while(std::getline(input_file, buffer)) {
		const std::vector<std::string> buffer_words = cafemol::library::split_String(buffer, ' ', '*');
	
		if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check " + input_name);
		if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
			is_reading_1stFrame = true;
			continue;
		}
		else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
			is_reading_1stFrame = false;
			continue;
		}
	
		if (is_reading_1stFrame) continue;
	//	else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && data_on_frame.empty()) continue;
		if ((buffer_words[0] == "Frame") && !is_reading_1stFrame) {
			all_data.row(i_frame) = data_on_frame;
			data_on_frame = Eigen::VectorXd::Zero(data_frame_x * data_frame_y);
			++i_frame;
			continue;
		}
	
	
		const std::size_t& data_coord_x = std::stoi(buffer_words[0]) - 1;
		const std::size_t& data_coord_y = std::stoi(buffer_words[1]) - 1;
	
		data_on_frame[data_coord_y * data_frame_x + data_coord_x] = 1.0;
	
	}
	
	all_data.row(i_frame) = data_on_frame;
//	if (!data_on_frame.empty()) {
//		all_data.row(i_frame) = data_on_frame;
	//	data_on_frame = Eigen::VectorXd::Zero(data_frame_x * data_frame_y);
//	}

	input_file.close();

	// ------------------------------------------------------------------------------------
	// read average

	sout("Read Average data");

	std::ifstream average_file(average_name, std::ios::in);

	buffer.clear();
	Eigen::VectorXd average_data(data_frame_x * data_frame_y);


	while (std::getline(average_file, buffer)) {
		const std::vector<std::string> buffer_words = cafemol::library::split_String(buffer, ' ', '*');
		const std::size_t& data_coord_x = std::stoi(buffer_words[0]) - 1;
		const std::size_t& data_coord_y = std::stoi(buffer_words[1]) - 1;
		const double& data_coord_z = std::stod(buffer_words[2]) - 1;

		average_data[data_coord_y * data_frame_x + data_coord_x] = data_coord_z;
	}
	average_file.close();

	// ------------------------------------------------------------------------------------
	// calc cross covariance matrix

	sout("Sum Cross Covariance Matrix");

	Eigen::MatrixXd&& sum_cross_covariance_matrix = Eigen::MatrixXd::Zero(all_data.cols(), all_data.cols());

	for (std::size_t i_row = 0; i_row < all_data.rows(); ++i_row) {
		sum_cross_covariance_matrix += (all_data.row(i_row).transpose() - average_data) * (all_data.row(i_row) - average_data.transpose());
	}

	// -------------------------------------------------------------------------------------
	// output
	
	sout("Output to " + output_name);
	std::ofstream ofs(output_name, std::ios::out | std::ios::binary);

	std::int32_t block_size = sizeof(int) * 2;
	int row_size = all_data.cols();
	int col_size = all_data.cols();
	ofs.write(reinterpret_cast<char*>(&block_size), sizeof(std::int32_t));
	ofs.write(reinterpret_cast<char*>(&row_size), sizeof(int));
	ofs.write(reinterpret_cast<char*>(&col_size), sizeof(int));
	ofs.write(reinterpret_cast<char*>(&block_size), sizeof(std::int32_t));

	std::int32_t mat_data_size = sizeof(double) * sum_cross_covariance_matrix.size();

	ofs.write(reinterpret_cast<char*>(&mat_data_size), sizeof(std::int32_t));

	for (int i_mat_datum = 0; i_mat_datum < sum_cross_covariance_matrix.size(); ++i_mat_datum) {
		ofs.write(reinterpret_cast<char*>(&sum_cross_covariance_matrix(i_mat_datum)), sizeof(double));
	}

	ofs.write(reinterpret_cast<char*>(&mat_data_size), sizeof(std::int32_t));

	ofs.close();

	return 0;
}


