#include<StandardOutput.hpp>
#include<ErrorMessage.hpp>
#include<OtherFormat/ContactStateReader.hpp>
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<memory>
#include<vector>
#include<array>
#include<sstream>
#include<iomanip>


int main(int argc, char* argv[]) {

	const int max_argc1 = 9;
	const int max_argc2 = 10;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if ((argc != max_argc1) && (argc != max_argc2)) eout("too much or less arguments.");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();
	bool is_set_frame_size;
	if (argc == max_argc1) is_set_frame_size = false;
	else if (argc == max_argc2) is_set_frame_size = true;

	const std::string& input_prefix = argv[1];
	const std::string& s_file_num = argv[2];
	const std::string& input_suffix = argv[3];
	const std::string& output_name = argv[4];
	const std::size_t& x_min_lim = std::stoi(argv[5]);
	const std::size_t& x_max_lim = std::stoi(argv[6]);
	const std::size_t& y_min_lim = std::stoi(argv[7]);
	const std::size_t& y_max_lim = std::stoi(argv[8]);

	const std::size_t& file_num = std::stoi(s_file_num);
	const std::size_t& file_num_digit = s_file_num.size();
	const std::array<std::size_t, 2> count_range_x = {x_min_lim, x_max_lim};
	const std::array<std::size_t, 2> count_range_y = {y_min_lim, y_max_lim};
	const std::size_t x_range = count_range_x[1] - count_range_x[0] + 1;
	const std::size_t y_range = count_range_y[1] - count_range_y[0] + 1;

	const std::size_t& BLOCK_SIZE = 90;

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Read Contact State Data", "");

	std::ifstream skip_file("SkipNum", std::ios::in);
	std::vector<std::size_t> skip_number;
	if (skip_file.is_open()) {
		sout("SkipNum has read.");
		std::string buffer;
		while (std::getline(skip_file, buffer)) {
			skip_number.push_back(std::stoi(buffer));
		}
	}
	skip_file.close();

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();

	std::vector<std::vector<std::array<int, 2>>> all_frames;
	std::vector<std::vector<std::size_t>> count_map(x_range);
	for (std::size_t idx = 0; idx < x_range; ++idx) {
		count_map[idx].resize(y_range, 0);
	}

	for (std::size_t file_index = 1; file_index <= file_num; ++file_index) {
		bool is_skipped = false;
		for (std::size_t skip_idx = 0; skip_idx < skip_number.size(); ++skip_idx) {
			is_skipped = (is_skipped || (skip_number[skip_idx] == file_index));
		}
		if (is_skipped) {
			sout("File Index " + std::to_string(file_index) + " skipped", "");
			continue;
		}

		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;

		std::string input_name(buffer_filename.str());
		contact_state_reader->read_ContactStatesWOIni(all_frames, input_name);
	}

	sout("calculating contacts");
	for (const std::vector<std::array<int, 2>>& frame : all_frames) {
		for (const std::array<int, 2>& contact_coordinates : frame) {
			const std::size_t& coords_x = static_cast<std::size_t>(contact_coordinates[0]);
			const std::size_t& coords_y = static_cast<std::size_t>(contact_coordinates[1]);
			if (((coords_x < x_min_lim) || (coords_x > x_max_lim)) || 
				((coords_y < y_min_lim) || (coords_y > y_max_lim))) continue;

			count_map[coords_x - x_min_lim][coords_y - y_min_lim] += 1;
		}
	}
	double frame_size;
	if (!is_set_frame_size) frame_size = static_cast<double>(all_frames.size());
	else {
		frame_size = std::stod(argv[9]);
	}

	sout("", "Output the results to " + output_name);
	std::ofstream ofs(output_name, std::ios::out);
	
	for (std::size_t x_idx = 0; x_idx < x_range; ++x_idx) {
		for (std::size_t y_idx = 0; y_idx < y_range; ++y_idx) {
			ofs << x_idx + x_min_lim << " "
				<< y_idx + y_min_lim << " "
				<< static_cast<double>(count_map[x_idx][y_idx]) / frame_size << std::endl;
		}
	}

	ofs.close();


	return 0;
}
