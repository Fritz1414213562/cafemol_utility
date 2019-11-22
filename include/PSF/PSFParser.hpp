#ifndef PSF_PARSER_HPP
#define PSF_PARSER_HPP
#include<fstream>
#include<iostream>
#include<array>
#include<vector>


namespace cafemol {

class PSFParser {

public:
	PSFParser() = default;
	PSFParser(const std::string& filename);
	~PSFParser() = default;
	void load_file();
	void close_file() {if (ifs.is_open()) ifs.close();}


protected:

	std::ifstream ifs;

	void load_file(const std::string& filename);

	int get_BlockNum(const std::string& line);
	std::string get_BlockKind(const std::string& line);

	int get_AtomID(const std::string& line);
	std::string get_SegmentName(const std::string& line);
	int get_ResidueNum(const std::string& line);
	std::string get_ResidueName(const std::string& line);
	std::string get_AtomName(const std::string& line);
	std::string get_AtomType(const std::string& line);
	float get_Charge(const std::string& line);
	float get_Mass(const std::string& line);

private:
	
	std::string input_name;

	enum COLUMN_MEANINGS {
		// header
		BLOCK_NUMBER_BEGIN = 0,
		BLOCK_NUMBER_END = 8,
		BLOCK_KIND_BEGIN = 9,
		BLOCK_KIND_END = 15,

		// contents in ATOM block
		ATOM_ID_BEGIN = 0,
		ATOM_ID_END = 8,
		SEGMENT_NAME_BEGIN = 9,
		SEGMENT_NAME_END = 13,
		RESIDUE_ID_BEGIN = 14,
		RESIDUE_ID_END = 18,
		RESIDUE_NAME_BEGIN = 19,
		RESIDUE_NAME_END = 23,
		ATOM_NAME_BEGIN = 24,
		ATOM_NAME_END = 27,
		ATOM_TYPE_BEGIN = 28,
		ATOM_TYPE_END = 31,
		CHARGE_BEGIN = 34,
		CHARGE_END = 44,
		MASS_BEGIN = 48,
		MASS_END = 58,
	};

	template<COLUMN_MEANINGS start, COLUMN_MEANINGS end>
	std::string readline_rangeof(const std::string& line) {
		int line_length = end - start;
		std::string res = line.substr(start, line_length);
		return res;
	}

	std::string trim_Space(const std::string& String);

};
}

#endif /* PSF_PARSER_HPP */
