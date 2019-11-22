#include"PSFParser.hpp"


cafemol::PSFParser::PSFParser(const std::string& filename) {
	input_name = filename;
	ifs.open(filename);
}

void cafemol::PSFParser::load_file(const std::string& filename) {
	close_file();
	input_name = filename;
	ifs.open(filename);
}


void cafemol::PSFParser::load_file() {
	close_file();
	ifs.open(input_name);
}


std::string cafemol::PSFParser::trim_Space(const std::string& String) {
	std::string res = "";
	for (const char& Char : String) {
		if (Char == ' ') continue;
		res += Char;
	}
	return res;
}


//std::string cafemol::PSFParser::readline_rangeof(const std::string& line, const int begin, const int end) {

//	std::string res;
//	int string_size = end - begin;
//	res = line.substr(begin, string_size);
//	return res;
//}


int cafemol::PSFParser::get_BlockNum(const std::string& line) {
	std::string buf = readline_rangeof<BLOCK_NUMBER_BEGIN, BLOCK_NUMBER_END>(line);
	int res = stoi(buf);
	return res;
}


std::string cafemol::PSFParser::get_BlockKind(const std::string& line) {
	std::string buf = readline_rangeof<BLOCK_KIND_BEGIN, BLOCK_KIND_END>(line);
	std::string res = trim_Space(buf);
	return res;
}


int cafemol::PSFParser::get_AtomID(const std::string& line) {
	std::string buf = readline_rangeof<ATOM_ID_BEGIN, ATOM_ID_END>(line);
	int res = stoi(buf);
	return res;
}


std::string cafemol::PSFParser::get_SegmentName(const std::string& line) {
	std::string buf = readline_rangeof<SEGMENT_NAME_BEGIN, SEGMENT_NAME_END>(line);
	std::string res = trim_Space(buf);
	return res;
}


int cafemol::PSFParser::get_ResidueNum(const std::string& line) {
	std::string buf = readline_rangeof<RESIDUE_ID_BEGIN, RESIDUE_ID_END>(line);
	int res = stoi(buf);
	return res;
}


std::string cafemol::PSFParser::get_ResidueName(const std::string& line) {
	std::string buf = readline_rangeof<RESIDUE_NAME_BEGIN, RESIDUE_NAME_END>(line);
	std::string res = trim_Space(buf);
	return res;
}


std::string cafemol::PSFParser::get_AtomName(const std::string& line) {
	std::string buf = readline_rangeof<ATOM_NAME_BEGIN, ATOM_NAME_END>(line);
	std::string res = trim_Space(buf);
	return res;
}


std::string cafemol::PSFParser::get_AtomType(const std::string& line) {
	std::string buf = readline_rangeof<ATOM_TYPE_BEGIN, ATOM_TYPE_END>(line);
	std::string res = trim_Space(buf);
	return res;
}


float cafemol::PSFParser::get_Charge(const std::string& line) {
	std::string buf = readline_rangeof<CHARGE_BEGIN, CHARGE_END>(line);
	float res = stof(buf);
	return res;
}


float cafemol::PSFParser::get_Mass(const std::string& line) {
	std::string buf = readline_rangeof<MASS_BEGIN, MASS_END>(line);
	float res = stof(buf);
	return res;
}
