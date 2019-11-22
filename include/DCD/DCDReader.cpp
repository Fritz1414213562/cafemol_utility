#include"DCDReader.hpp"

cafemol::DCDReader::DCDReader(const std::string input_file_name) : DCDParser(input_file_name) {
	frame_num = 0;
	atom_num = 0;
}


std::array<std::vector<float>, 3> cafemol::DCDReader::get_SnapShot_at(const int frame_in_snapshot) {
	
	// read header, and store frame_num and atom_num
	read_num_frame_and_atom();
	std::array<std::vector<float>, 3> result;

	if (frame_in_snapshot > frame_num) {
		std::cout << "The frame '" << frame_in_snapshot << "' is out of range." << std::endl;
		std::exit(1);
	}

	for (std::size_t iframe = 0; iframe < frame_in_snapshot; ++iframe) {
		result = read_xyz(atom_num);
	}

	// reload file.
	close_file();
	open_file();

	return result;
}
