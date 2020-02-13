#include<OtherFormat/ContactStateReader.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<sstream>
#include<iomanip>
#include<fstream>
#include<memory>
#include<array>
#include<vector>
#include<string>


int main(int argc, char *argv[]) {

	const int& max_argc = 8;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_prefix = argv[1];
	const std::size_t& file_num = std::stoi(argv[2]);
	const std::string& input_suffix = argv[3];
	const std::string& state_definition_prefix = argv[4];
	const std::string& state_definition_suffix = argv[5];
	const std::string& output_prefix = argv[6];
	const std::string& output_suffix = argv[7];

	const std::size_t& BLOCK_SIZE = 90;

	const std::string& s_file_num = argv[2];
	const std::size_t& file_num_digit = s_file_num.size();

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Reading Contact State Data");

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

	for (std::size_t file_index = 1; file_index <= file_num; ++file_index) {
		bool is_skipped = false;
		for (std::size_t skip_idx = 0; skip_idx < skip_number.size(); ++skip_idx) {
			is_skipped = (is_skipped || (skip_number[skip_idx] == file_index));
		}
		if (is_skipped) {
			sout("File Index " + std::to_string(file_index) + " skipped.", "");
			continue;
		}

		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;

		std::stringstream buffer_statedef;
		buffer_statedef << state_definition_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< state_definition_suffix;

		std::string input_name(buffer_filename.str());
		std::string state_def_name(buffer_statedef.str());

		std::ifstream state_file(state_def_name, std::ios::in);
		std::string buffer;
		std::vector<std::size_t> state_frames;

		while (std::getline(state_file, buffer)) {
			std::size_t state_end_frame = std::stoi(buffer);
			state_frames.push_back(state_end_frame);
		}
		state_file.close();
		sout("The file, " + state_def_name + " has read.");


		const std::vector<std::vector<std::array<int, 2>>>& frames = contact_state_reader->read_ContactStatesWOIni(input_name);
		sout("The file, " + input_name + " has read.");
		sout("");

//		for (std::size_t idx = 0; idx < state_frames.size(); ++idx) {
//			if (state_frames[idx] == 0) {
//				state_frames[idx] = frames.size();
//			}
//		}
		state_frames.push_back(frames.size());

		sout("Contact States will be divided into " + std::to_string(state_frames.size()));
		for (const std::size_t& state_end_frame : state_frames) {
			sout(state_end_frame);
		}

		std::size_t divided_state_idx = 1;
		std::size_t frame_start_idx = 0;
		for (std::size_t state_end_frame : state_frames) {
			std::stringstream buffer_name;
			buffer_name << output_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< "_state"
						<< divided_state_idx
						<< output_suffix;

			const std::string& output_name(buffer_name.str());
			std::ofstream ofs(output_name, std::ios::out);

			for (; frame_start_idx < state_end_frame; ++frame_start_idx) {
				ofs << "Frame " << frame_start_idx + 1 << std::endl;
				for (const std::array<int, 2>& contact_coords : frames[frame_start_idx]) {
					ofs << contact_coords[0] << " " << contact_coords[1] << std::endl;
				}
			}
			ofs.close();
			++divided_state_idx;
		}



	}

	sout("Program finished");

	return 0;
}
