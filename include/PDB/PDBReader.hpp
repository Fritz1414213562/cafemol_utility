#pragma once
#include"PDBParser.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<array>


namespace cafemol {


class PDBReader : public PDBParser {

public:
	PDBReader() = default;
	~PDBReader() = default;
	PDBReader(const std::string& filename) : PDBParser(filename) {}
	void read_PDBData();

	int get_ResidueIDof(const int& atom_id);
	int get_AtomIDof(const int& local_residue_num, const std::string& chainID);
	int get_AtomNumof(const std::string& query);
	int get_AllAtomNum();
	int get_ResidueNumof(const std::string& chainID);

	std::array<float, 3> get_Coordinateof(const int& atom_id);
	std::vector<std::array<float, 3>> get_Coordinateof(const int& atom_id_begin, const int& atom_id_end);

private:
	std::vector<int> residue_id_container;
	std::vector<int> atom_id_container;
	std::vector<int> chain_begin_indices;
	std::vector<std::array<float, 3>> coordinates_container;
	std::string chainIDs;


	void combine_AtomID2ResiID();
	void read_AllChainID();
	void read_Coordinates();

	bool is_ATOMColumn(const std::string& line) {
		std::string rec_type = read_RecordName(line);
		bool res = ((rec_type == "ATOM  ") || (rec_type == "TER"));
		return res;
	}

	void is_OpenInputStream() {
		if (!input_file.is_open()) {
			std::cerr << "Error: The file has not been yet opened." << std::endl;
			std::exit(1);
		}
	}

	int getIndexof(const int& atom_id) {
		int res = 0;
		for (const int& id : atom_id_container) {
			if (id == atom_id) break;
			else ++res;
		}

		// when atom id does not exist in container.
		int index_out_of_container = atom_id_container.size();
		if (res >= index_out_of_container) {
			std::cerr << "Error: The Atom ID '" << atom_id << "' does not exist." << std::endl;
			std::exit(1);
		}
		return res;
	}

};
}
