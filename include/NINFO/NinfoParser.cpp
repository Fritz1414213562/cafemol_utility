#include"NinfoParser.hpp"


cafemol::NinfoParser::NinfoParser(const std::string& file_name) {
	input_name = file_name;
	input_file.open(file_name, std::ios::in);
}


void cafemol::NinfoParser::load_file(const std::string& file_name) {
	close_file();
	input_name = file_name;
	input_file.open(file_name, std::ios::in);
}


void cafemol::NinfoParser::load_file() {
	close_file();
	input_file.open(input_name, std::ios::in);
}


void cafemol::NinfoParser::close_file() {if (input_file.is_open()) input_file.close();}


std::vector<std::string> cafemol::NinfoParser::split_Line(const std::string& line, const char& delimiter) {

	std::vector<std::string> result;
	std::string buf;
	for (const char& character : line) {

		if ((buf.empty()) && (character == delimiter)) continue;
		else if ((!buf.empty()) && (character == delimiter)) {
			result.push_back(buf);
			buf.clear();
			continue;
		}

		buf.push_back(character);
	}
	if (!buf.empty()) result.push_back(buf);

	return result;
}


std::string cafemol::NinfoParser::get_Headof(const std::string& line) {
	std::vector<std::string> line_vec = split_Line(line, ninfo_delimiter);
	std::string res;
	if (line_vec.empty()) return res;
	else if (!line_vec.empty())  return line_vec[0];
}


bool cafemol::NinfoParser::is_Block(const std::string& line) {

	bool result;
	std::vector<std::string> line_vec = split_Line(line, ninfo_delimiter);
	if (line_vec[0] == "<<<<") result = true;
	else result = false;

	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::bond> cafemol::NinfoParser::read_Line(const std::string& line) {
	std::vector<std::string> buf = split_Line(line, ninfo_delimiter);
	ninfo_data_type::nativeinfo_tuple<enum_types::bond> result;

	result.line_kind_name = buf[0];
	std::get<0>(result.line_data) = std::stoi(buf[1]);
	std::get<1>(result.line_data) = {std::stoi(buf[2]), std::stoi(buf[3])};
	std::get<2>(result.line_data) = {std::stoi(buf[4]), std::stoi(buf[5]), 
									 std::stoi(buf[6]),	std::stoi(buf[7])};
	std::get<3>(result.line_data) = {std::stof(buf[8]), std::stof(buf[9]),
									 std::stof(buf[10]), std::stof(buf[11])};
	
	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::angl> cafemol::NinfoParser::read_Line(const std::string& line) {
	std::vector<std::string> buf = split_Line(line, ninfo_delimiter);
	ninfo_data_type::nativeinfo_tuple<enum_types::angl> result;

	result.line_kind_name = buf[0];
	std::get<0>(result.line_data) = std::stoi(buf[1]);
	std::get<1>(result.line_data) = {std::stoi(buf[2]), std::stoi(buf[3])};
	std::get<2>(result.line_data) = {std::stoi(buf[4]), std::stoi(buf[5]), 
									 std::stoi(buf[6]),	std::stoi(buf[7]),
									 std::stoi(buf[8]),	std::stoi(buf[9])};
	std::get<3>(result.line_data) = {std::stof(buf[10]), std::stof(buf[11]),
									 std::stof(buf[12]), std::stof(buf[13])};
	
	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::aicg13> cafemol::NinfoParser::read_Line(const std::string& line) {
	std::vector<std::string> buf = split_Line(line, ninfo_delimiter);
	ninfo_data_type::nativeinfo_tuple<enum_types::aicg13> result;

	result.line_kind_name = buf[0];
	std::get<0>(result.line_data) = std::stoi(buf[1]);
	std::get<1>(result.line_data) = {std::stoi(buf[2]), std::stoi(buf[3])};
	std::get<2>(result.line_data) = {std::stoi(buf[4]), std::stoi(buf[5]), 
									 std::stoi(buf[6]),	std::stoi(buf[7]),
									 std::stoi(buf[8]),	std::stoi(buf[9])};
	std::get<3>(result.line_data) = {std::stof(buf[10]), std::stof(buf[11]),
									 std::stof(buf[12]), std::stof(buf[13]),
									 std::stof(buf[14])};
	
	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::dihd> cafemol::NinfoParser::read_Line(const std::string& line) {
	std::vector<std::string> buf = split_Line(line, ninfo_delimiter);
	ninfo_data_type::nativeinfo_tuple<enum_types::dihd> result;

	result.line_kind_name = buf[0];
	std::get<0>(result.line_data) = std::stoi(buf[1]);
	std::get<1>(result.line_data) = {std::stoi(buf[2]), std::stoi(buf[3])};
	std::get<2>(result.line_data) = {std::stoi(buf[4]), std::stoi(buf[5]), 
									 std::stoi(buf[6]),	std::stoi(buf[7]),
									 std::stoi(buf[8]),	std::stoi(buf[9]),
									 std::stoi(buf[10]), std::stoi(buf[11])};
	std::get<3>(result.line_data) = {std::stof(buf[12]), std::stof(buf[13]),
									 std::stof(buf[14]), std::stof(buf[15]),
									 std::stof(buf[16])};
	
	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::aicgdih> cafemol::NinfoParser::read_Line(const std::string& line) {
	std::vector<std::string> buf = split_Line(line, ninfo_delimiter);
	ninfo_data_type::nativeinfo_tuple<enum_types::aicgdih> result;

	result.line_kind_name = buf[0];
	std::get<0>(result.line_data) = std::stoi(buf[1]);
	std::get<1>(result.line_data) = {std::stoi(buf[2]), std::stoi(buf[3])};
	std::get<2>(result.line_data) = {std::stoi(buf[4]), std::stoi(buf[5]), 
									 std::stoi(buf[6]),	std::stoi(buf[7]),
									 std::stoi(buf[8]),	std::stoi(buf[9]),
									 std::stoi(buf[10]), std::stoi(buf[11])};
	std::get<3>(result.line_data) = {std::stof(buf[12]), std::stof(buf[13]),
									 std::stof(buf[14]), std::stof(buf[15]),
									 std::stof(buf[16])};
	
	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::contact> cafemol::NinfoParser::read_Line(const std::string& line) {
	std::vector<std::string> buf = split_Line(line, ninfo_delimiter);
	ninfo_data_type::nativeinfo_tuple<enum_types::contact> result;

	result.line_kind_name = buf[0];
	std::get<0>(result.line_data) = std::stoi(buf[1]);
	std::get<1>(result.line_data) = {std::stoi(buf[2]), std::stoi(buf[3])};
	std::get<2>(result.line_data) = {std::stoi(buf[4]), std::stoi(buf[5]), 
									 std::stoi(buf[6]),	std::stoi(buf[7])};
	std::get<3>(result.line_data) = {std::stof(buf[8]), std::stof(buf[9]),
									 std::stof(buf[10]), std::stof(buf[11])};
	
	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::pdns> cafemol::NinfoParser::read_Line(const std::string& line) {
	std::vector<std::string> buf = split_Line(line, ninfo_delimiter);
	ninfo_data_type::nativeinfo_tuple<enum_types::pdns> result;

	result.line_kind_name = buf[0];
	std::get<0>(result.line_data) = std::stoi(buf[1]);
	std::get<1>(result.line_data) = {std::stoi(buf[2])};
	std::get<2>(result.line_data) = {std::stoi(buf[3]), std::stoi(buf[4])};
	std::get<3>(result.line_data) = {std::stof(buf[5]), std::stof(buf[6]),
									 std::stof(buf[7]), std::stof(buf[8])};
	
	return result;
}


template<>
cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::other> cafemol::NinfoParser::read_Line(const std::string& line) {
	ninfo_data_type::nativeinfo_tuple<enum_types::other> result;

	result.line = line;
	return result;
}


