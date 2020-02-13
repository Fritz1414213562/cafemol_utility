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


int main(int argc, char* argv[]) {

	const int max_argc = 7;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments.");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_name = argv[1];
	const std::string& output_name = argv[2];
	const std::size_t& x_min_lim = std::stoi(argv[3]);
	const std::size_t& x_max_lim = std::stoi(argv[4]);
	const std::size_t& y_min_lim = std::stoi(argv[5]);
	const std::size_t& y_max_lim = std::stoi(argv[6]);
	const std::array<std::size_t, 2> count_range_x = {x_min_lim, x_max_lim};
	const std::array<std::size_t, 2> count_range_y = {y_min_lim, y_max_lim};
	const std::size_t x_range = count_range_x[1] - count_range_x[0] + 1;
	const std::size_t y_range = count_range_y[1] - count_range_y[0] + 1;

	const std::size_t& BLOCK_SIZE = 90;


	std::vector<std::vector<std::array<int, 2>>> all_frames;
	std::vector<std::vector<std::size_t>> count_map(x_range);
	for (std::size_t idx = 0; idx < x_range; ++idx) {
		count_map[idx].resize(y_range, 0);
	}

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Read Contact State Data", "");
	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();
	contact_state_reader->read_ContactStatesWOIni(all_frames, input_name);

	sout("Count contacts");
	for (const std::vector<std::array<int, 2>>& frame : all_frames) {
		for (const std::array<int, 2>& contact_coordinates : frame) {
			const std::size_t& coords_x = static_cast<std::size_t>(contact_coordinates[0]);
			const std::size_t& coords_y = static_cast<std::size_t>(contact_coordinates[1]);
			if (((coords_x < x_min_lim) || (coords_x > x_max_lim)) || 
				((coords_y < y_min_lim) || (coords_y > y_max_lim))) continue;

			count_map[coords_x - x_min_lim][coords_y - y_min_lim] += 1;
		}
	}

	sout("", "Output the results to " + output_name);
	std::ofstream ofs(output_name, std::ios::out);
	
	for (std::size_t x_idx = 0; x_idx < x_range; ++x_idx) {
		for (std::size_t y_idx = 0; y_idx < y_range; ++y_idx) {
			ofs << x_idx + x_min_lim << " "
				<< y_idx + y_min_lim << " "
				<< count_map[x_idx][y_idx] << std::endl;
		}
	}

	ofs.close();


	return 0;
}
