#include<OtherFormat/MatrixFileReader.hpp>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<UtilFunc.hpp>
#include<vector>
#include<string>
#include<memory>
#include<fstream>
#include<sstream>
#include<iomanip>


int main(int argc, char* argv[]) {

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != 8) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_name = argv[1];
	const std::string& axes_name = argv[2];
	const std::string& output_name = argv[3];
	const std::size_t& max_component_num = std::stoi(argv[4]);
	const std::size_t& data_frame_x = std::stoi(argv[5]);
	const std::size_t& data_frame_y = std::stoi(argv[6]);
	const std::size_t& frame_size = std::stoi(argv[7]);

	const std::size_t& data_size = data_frame_x * data_frame_y;
	if (data_size < max_component_num) eout("Max Componenet number is larger than the original data size");

	std::string buffer;
	bool is_reading_1stFrame = false;
	Eigen::MatrixXd all_data = Eigen::MatrixXd::Zero(frame_size, data_size);
	Eigen::VectorXd data_on_frame = Eigen::VectorXd::Zero(data_size);

	sout("Read Contact State Data");
	std::ifstream input_file(input_name, std::ios::in);
	std::size_t i_frame = 0;

	while(std::getline(input_file, buffer)) {
		const std::vector<std::string> buffer_words = cafemol::library::split_String(buffer, ' ', '*');

		if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", input_name);
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
			data_on_frame = Eigen::VectorXd::Zero(data_size);
			++i_frame;
			continue;
		}


		const std::size_t& data_coord_x = std::stoi(buffer_words[0]) - 1;
		const std::size_t& data_coord_y = std::stoi(buffer_words[1]) - 1;

		data_on_frame[data_coord_y * data_frame_x + data_coord_x] = 1.0;

	}

	if (i_frame < frame_size) {
		all_data.row(i_frame) = data_on_frame;
	}
	

	input_file.close();
	sout(".... Done");

	sout("Read Principal Componenet Axes Data");
	std::unique_ptr<cafemol::MatrixFileReader> matrix_reader = std::make_unique<cafemol::MatrixFileReader>();

	const Eigen::MatrixXd& eigen_values = matrix_reader->read_Matrix(axes_name);
	if (eigen_values.rows() != data_size) eout("The dimension of pc axes is not consistent with input data size");

	Eigen::MatrixXd pc_axes(data_size, max_component_num);
	sout("PC size = " + std::to_string(max_component_num));
	
	for (std::size_t i_axis = 0; i_axis < max_component_num; ++i_axis) {
		pc_axes.col(i_axis) = eigen_values.col((eigen_values.cols() - 1) - i_axis);
	}
	sout(".... Done");

	sout("Projection data to pc axes");

	const Eigen::MatrixXd& result = all_data * pc_axes;
	sout("Output to " + output_name);

	std::ofstream ofs(output_name, std::ios::out);

	ofs << result << std::endl;
	ofs.close();

	return 0;
}
