#include"PDBParser.hpp"


	// parse method
std::string cafemol::PDBParser::read_RecordName(const std::string& line) {
	std::string buf = readline_range_of(line, RECORD_NAME_BEGIN, TER_END);
	if (buf == "TER") return "TER";
	else return readline_range_of(line, RECORD_NAME_BEGIN, RECORD_NAME_END);
}


int cafemol::PDBParser::read_SerialNum(const std::string& line) {
	std::string buf = readline_range_of(line, SERIAL_NUM_BEGIN, SERIAL_NUM_END);
	int res = stoi(buf);
	return res;
}


std::string cafemol::PDBParser::read_AtomName(const std::string& line) {
	return readline_range_of(line, ATOM_NAME_BEGIN, ATOM_NAME_END);
}


char cafemol::PDBParser::read_AltLoc(const std::string& line) {
	std::string buf = readline_range_of(line, ALT_LOC, ALT_LOC);
	char res = *buf.c_str();
	return res;
}


std::string cafemol::PDBParser::read_ResidueName(const std::string& line) {
	return readline_range_of(line, RESIDUE_NAME_BEGIN, RESIDUE_NAME_END);
}


std::string cafemol::PDBParser::read_ChainID(const std::string& line) {
	return readline_range_of(line, CHAIN_ID, CHAIN_ID);
}


int cafemol::PDBParser::read_ResiSerialNum(const std::string& line) {
	std::string buf = readline_range_of(line, RESI_SEQ_NUM_BEGIN, RESI_SEQ_NUM_END);
	int res = stoi(buf);
	return res;
}


char cafemol::PDBParser::read_InsertResiCode(const std::string& line) {
	std::string buf = readline_range_of(line, INSERTI_RESI_NUM, INSERTI_RESI_NUM);
	char res = *buf.c_str();
	return res;
}


float cafemol::PDBParser::read_CoordinateX(const std::string& line) {
	std::string buf = readline_range_of(line, CORD_X_BEGIN, CORD_X_END);
	float res = stof(buf);
	return res;
}


float cafemol::PDBParser::read_CoordinateY(const std::string& line) {
	std::string buf = readline_range_of(line, CORD_Y_BEGIN, CORD_Y_END);
	float res = stof(buf);
	return res;
}


float cafemol::PDBParser::read_CoordinateZ(const std::string& line) {
	std::string buf = readline_range_of(line, CORD_Z_BEGIN, CORD_Z_END);
	float res = stof(buf);
	return res;
}


float cafemol::PDBParser::read_Occupancy(const std::string& line) {
	std::string buf = readline_range_of(line, OCCUPANCY_BEGIN, OCCUPANCY_END);
	float res = stof(buf);
	return res;
}


float cafemol::PDBParser::read_TempFactor(const std::string& line) {
	std::string buf = readline_range_of(line, TEMP_FACTOR_BEGIN, TEMP_FACTOR_END);
	float res = stof(buf);
	return res;
}


std::string cafemol::PDBParser::read_Element(const std::string& line) {
	return readline_range_of(line, ELEMENT_BEGIN, ELEMENT_END);
}


float cafemol::PDBParser::read_Charge(const std::string& line) {
	std::string buf = readline_range_of(line, CHARGE_BEGIN, CHARGE_END);
	float res = stof(buf);
	return res;
}

