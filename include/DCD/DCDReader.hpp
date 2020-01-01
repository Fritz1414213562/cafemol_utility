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
	std::array<std::vector<float>, 3> get_SnapShot_at(const int frame_in_snapshot);

};

}

#endif /* DCD_READER_HPP */
