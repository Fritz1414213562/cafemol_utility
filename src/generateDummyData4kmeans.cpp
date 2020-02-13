#include<OtherFormat/ContactStateReader.hpp>
#include<StandardOutput.hpp>
#include<ErrorMessage.hpp>
#include<random>
#include<string>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<memory>
#include<array>
#include<vector>


int main(int argc, char* argv[]) {

	const int& max_argc = 10;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments.");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string input_prefix = argv[1];
	const std::size_t file_num = std::stoi(argv[2]);
	const std::string input_suffix = argv[3];
	const std::string output_name = argv[4];

	const std::size_t random_seed = std::stoi(argv[5]);

	const std::size_t x_min = std::stoi(argv[6]);
	const std::size_t x_max = std::stoi(argv[7]);
	const std::size_t y_min = std::stoi(argv[8]);
	const std::size_t y_max = std::stoi(argv[9]);

	const std::size_t& BLOCK_SIZE = 90;

	const std::string& s_file_num = argv[2];
	const std::size_t& file_num_digit = s_file_num.size();

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();

	std::vector<std::size_t> all_frame_sizes;

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Reading ContactState Data");

	std::ifstream skip_file("SkipNum", std::ios::in);
	std::vector<std::size_t> skip_number;
	if (skip_file.is_open()) {
		sout("The file, SkipNum has read.");
		std::string buffer;
		while(std::getline(skip_file, buffer)) {
			skip_number.push_back(std::stoi(buffer));
		}
	}
	skip_file.close();



	for (std::size_t file_index = 1; file_index <= file_num; ++file_index) {
		bool is_skipped = false;
		for (std::size_t skip_index = 0; skip_index < skip_number.size(); ++skip_index) {
			is_skipped = (is_skipped || (skip_number[skip_index] == file_index));
		}
		if (is_skipped) {
			sout("File Index " + std::to_string(file_index) + " skipped.");
			continue;
		}

		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;
		std::string input_name(buffer_filename.str());
		contact_state_reader->read_ContactStatesFrameSizeWOIni(all_frame_sizes, input_name);
	}
	sout("", ".... Done");

	sout("Generating dummy data");

	std::mt19937_64 mt_engine_x(random_seed);
	std::mt19937_64 mt_engine_y(random_seed + 1);
	std::uniform_int_distribution<> dist_x(x_min, x_max);
	std::uniform_int_distribution<> dist_y(y_min, y_max);

	std::vector<std::vector<std::array<std::size_t, 2>>> result(all_frame_sizes.size());

	for (std::size_t i_frame = 0; i_frame < all_frame_sizes.size(); ++i_frame) {
		const std::size_t& frame_size = all_frame_sizes[i_frame];
		result[i_frame].resize(frame_size);
		for (std::size_t i_datum = 0; i_datum < frame_size; ++i_datum) {
			std::size_t datum_x = dist_x(mt_engine_x);
			std::size_t datum_y = dist_y(mt_engine_y);
			result[i_frame][i_datum] = {datum_x, datum_y};
		}
	}

	sout("Output to " + output_name);
	std::ofstream ofs(output_name, std::ios::out);

	std::size_t i_frame = 1;
	for (const std::vector<std::array<std::size_t, 2>>& dummy_data : result) {
		ofs << "Frame " << i_frame << std::endl;
		for (const std::array<std::size_t, 2>& dummy_datum : dummy_data) {
			ofs << dummy_datum[0] << " " << dummy_datum[1] << std::endl;
		}
		++i_frame;
	}

	sout("", ".... Done");
	sout.output_HyphenBlock("", BLOCK_SIZE);

	return 0;

}
