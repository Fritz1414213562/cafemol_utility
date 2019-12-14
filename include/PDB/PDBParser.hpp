// PDBParser
#ifndef PDB_PARSER_HPP
#define PDB_PARSER_HPP

#include<string>
#include<iostream>
#include<fstream>

// namespace
namespace cafemol {

class PDBParser {


public:

	PDBParser() = default;
	~PDBParser() = default;
	PDBParser(const std::string& filename) {load_file(filename);}

	void load_file(const std::string& file_name) {
		close_file();
		input_name = file_name;
		input_file.open(file_name, std::ios::in);
	}


protected:

	std::ifstream input_file;
	std::string input_name;

	// open / close file 
	void load_file() {
		close_file();
		input_file.open(input_name, std::ios::in);
	}

	void close_file() {if (input_file.is_open()) input_file.close();}

	// parse method
	std::string read_RecordName(const std::string& line);

	int read_SerialNum(const std::string& line);

	std::string read_AtomName(const std::string& line);

	char read_AltLoc(const std::string& line);

	std::string read_ResidueName(const std::string& line);

	std::string read_ChainID(const std::string& line);

	int read_ResiSerialNum(const std::string& line);

	char read_InsertResiCode(const std::string& line);

	float read_CoordinateX(const std::string& line);

	float read_CoordinateY(const std::string& line);

	float read_CoordinateZ(const std::string& line);

	float read_Occupancy(const std::string& line);

	float read_TempFactor(const std::string& line);

	std::string read_Element(const std::string& line);

	float read_Charge(const std::string& line);


private:

	// the number columns
	enum Column_numbers {
		RECORD_NAME_BEGIN = 0,
		TER_END = 2,
		RECORD_NAME_END = 5,
		SERIAL_NUM_BEGIN = 6,
		SERIAL_NUM_END = 10,
		ATOM_NAME_BEGIN = 12,
		ATOM_NAME_END = 15,
		ALT_LOC = 16,
		RESIDUE_NAME_BEGIN = 17,
		RESIDUE_NAME_END = 19,
		CHAIN_ID = 21,
		RESI_SEQ_NUM_BEGIN = 22,
		RESI_SEQ_NUM_END = 25,
		INSERTI_RESI_NUM = 26,
		CORD_X_BEGIN = 30,
		CORD_X_END = 37,
		CORD_Y_BEGIN = 38,
		CORD_Y_END = 45,
		CORD_Z_BEGIN = 46,
		CORD_Z_END = 53,
		OCCUPANCY_BEGIN = 54,
		OCCUPANCY_END = 59,
		TEMP_FACTOR_BEGIN = 60,
		TEMP_FACTOR_END = 65,
		ELEMENT_BEGIN = 76,
		ELEMENT_END = 77,
		CHARGE_BEGIN = 78,
		CHARGE_END = 79
	};

	std::string readline_range_of(const std::string& line, const int begin, const int end) {
		std::string res;
		int string_size = end - begin + 1;
		res = line.substr(begin, string_size);
		return res;
	}
};
}

#endif /* PDB_PARSER_HPP */
