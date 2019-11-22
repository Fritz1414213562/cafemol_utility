// This class can make pdb file from coordinates data.
#pragma once
#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<array>
#include<iomanip>

// namespace
namespace cafemol {

class PDBWriter {

public:
	
	PDBWriter() = default;
	~PDBWriter() = default;

	void makePDBFileFrom(const std::string& template_name, const std::string& out_name, const std::array<std::vector<float>, 3>& vecs);

private:
	
	void close_file(std::ifstream& input_file);

	// methods for making PDB file from a template PDB file.
	void replace_CodsInPDBLine(std::ofstream& ofs, const std::string& line, const std::array<float, 3>& vec);

	// read line from a template PDB file
	std::string read_LineBeforeCods(const std::string& line);
	std::string read_LineAfterCods(const std::string& line);
	int read_SerialNumof(const std::string& line);

	enum PDB_COLUMN_MEANING {
		// the column numbers of PDB format
		LINE_BEGIN = 0,
		TER_COLUMN = 2,
		REC_TYPE_ATOM_END = 5,
		SERIAL_NUM_BEGIN = 6,
		SERIAL_NUM_END = 10,
		BEFORE_COORDINATES = 29,
		AFTER_COORDINATES = 54,

		// the column size of coordinates
		COORDINATES_DECIMAL_SIZE = 3,
		COORDINATES_COLUMN_SIZE = 8,
	};

	std::string read_LineRangeof(const std::string& line, const int begin, const int end);
	
	bool is_AtomRow(const std::string& line);
	bool is_TERRow(const std::string& line);

	std::vector<std::array<float, 3>> transpose_Mat(const std::array<std::vector<float>, 3>& matrix);
};


// -----------------------------------------------------------------------------------------------
// inline function

inline std::string cafemol::PDBWriter::read_LineRangeof(const std::string& line, const int begin, const int end) {
		std::string res;
		int string_size = end - begin + 1;
		res = line.substr(begin, string_size);
		return res;
	}
}


inline bool cafemol::PDBWriter::is_AtomRow(const std::string& line) {
	std::string rec_type = read_LineRangeof(line, LINE_BEGIN, REC_TYPE_ATOM_END);
	bool res = (rec_type == "ATOM  ");
	return res;
}


inline bool cafemol::PDBWriter::is_TERRow(const std::string& line) {
	std::string rec_type = read_LineRangeof(line, LINE_BEGIN, TER_COLUMN);
	bool res = (rec_type == "TER");
	return res;
}
