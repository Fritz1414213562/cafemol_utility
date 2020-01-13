// This class can read the coordinates in a snapshot.
#ifndef DCD_READER_HPP
#define DCD_READER_HPP

#include"DCDParser.hpp"
#include<iostream>
#include<fstream>
#include<vector>
#include<array>

namespace cafemol {

class DCDReader : public DCDParser {

public :
	DCDReader(const std::string& input_file_name);

	std::array<std::vector<float>, 3> get_SnapShot_at(const std::size_t& frame_in_snapshot);

	std::vector<std::array<std::vector<float>, 3>> get_Trajectory();

};

}

#endif /* DCD_READER_HPP */
