#include"PDBReader.hpp"


void cafemol::PDBReader::read_PDBData() {
	std::cout << "----------------------------------------------------------------------------------" << std::endl;
	combine_AtomID2ResiID();
	std::cout << "Read Atom ID and Residue ID" << std::endl;
	read_AllChainID();
	std::cout << "Read chain ID" << std::endl;
	read_Coordinates();
	std::cout << "Read the coordinates" << std::endl;
	std::cout << "----------------------------------------------------------------------------------" << std::endl;
}


void cafemol::PDBReader::combine_AtomID2ResiID() {

	is_OpenInputStream();

	atom_id_container.clear();
	residue_id_container.clear();
	chain_begin_indices.clear();
	chain_begin_indices.push_back(0);
	// initialize line
	std::string line;

	int chain_begin_atom_index = 0;

	while (getline(input_file, line)) {
		if (!is_ATOMColumn(line)) continue;

		std::string record_type = read_RecordName(line);
		if (record_type == "TER") {
			chain_begin_indices.push_back(chain_begin_atom_index);
			continue;
		}

		int atom_id = read_SerialNum(line);
		int residue_id = read_ResiSerialNum(line);

		atom_id_container.push_back(atom_id);
		residue_id_container.push_back(residue_id);
		++chain_begin_atom_index;
	}

	// reload file stream for the next use
	load_file();
}


void cafemol::PDBReader::read_AllChainID() {
	is_OpenInputStream();

	chainIDs = "";
	// initialize line
	std::string line;
	
	std::string chain_id;

	while (getline(input_file, line)) {
		if (!is_ATOMColumn(line)) continue;

		std::string record_type = read_RecordName(line);
		if (record_type == "TER") {
			chainIDs += chain_id;
			chain_id.clear();
			continue;
		}

		if (chain_id.empty()) {
			chain_id = read_ChainID(line);
			continue;
		}

		else if (chain_id != read_ChainID(line)) {
			std::cerr << "Error: This PDB file doesn't have 'TER' column between different chains" << std::endl;
			std::exit(1);
		}
	}

	// reload file stream for the next use
	load_file();
}


void cafemol::PDBReader::read_Coordinates() {
	is_OpenInputStream();

	coordinates_container.clear();

	std::string line;
	while (std::getline(input_file, line)) {
		if (!is_ATOMColumn(line)) continue;

		std::string record_type = read_RecordName(line);
		if (record_type == "TER") continue;

		float coordinate_x = read_CoordinateX(line);
		float coordinate_y = read_CoordinateY(line);
		float coordinate_z = read_CoordinateZ(line);
		coordinates_container.push_back({coordinate_x, coordinate_y, coordinate_z});
	}
	// reload file stream for the next use
	load_file();
}


int cafemol::PDBReader::get_ResidueIDof(const int& atom_id) {

	// Judge whether the containers have their information.
	if (residue_id_container.empty() || atom_id_container.empty()) {
		std::cerr << "Error: The container of resi or atom is empty." << std::endl;
		std::exit(1);
	}

	// search atom_id in atom_id_container.
	int atom_id_index = getIndexof(atom_id);

	int result = residue_id_container[atom_id_index];
	return result;
}


std::array<float, 3> cafemol::PDBReader::get_Coordinateof(const int& atom_id) {

	// judge whether the containers have their information.
	if (coordinates_container.empty() || atom_id_container.empty()) {
		std::cerr << "Error: The container of coordinates or atom IDs is empty." << std::endl;
		std::exit(1);
	}

	// search atom_id in atom_id_container.
	int atom_id_index = getIndexof(atom_id);
	std::array<float, 3> result = coordinates_container[atom_id_index];
	return result;
}


std::vector<std::array<float, 3>> cafemol::PDBReader::get_Coordinateof(const int& atom_id_begin, const int& atom_id_end) {

	// judge whether the containers have their information.
	if (coordinates_container.empty() || atom_id_container.empty()) {
		std::cerr << "Error: The container of coordinates or atom IDs is empty." << std::endl;
		std::exit(1);
	}

	// search atom_id in atom_id_container.
	int atom_id_begin_index = getIndexof(atom_id_begin);
	int atom_id_end_index = getIndexof(atom_id_end);

	std::vector<std::array<float, 3>> result;
	for (std::size_t idx = atom_id_begin_index; idx <= atom_id_end_index; ++idx) {
		result.push_back(coordinates_container[idx]);
	}
	return result;
}


int cafemol::PDBReader::get_AtomIDof(const int& local_residue_num, const std::string& chainID) {
	
	// judge whether the string have thein information
	if (chainIDs.empty()) {
		std::cerr << "Error: The string of Chain ID is empty." << std::endl;
		std::exit(1);
	}

	// search 
	int chain_id_index = chainIDs.find(chainID);
	int chain_index_begin = chain_begin_indices[chain_id_index];
	int chain_index_end = chain_begin_indices.back();
	if (chain_id_index < chainIDs.size() - 1) {
		chain_index_end = chain_begin_indices[chain_id_index + 1];
	}

	int idx_of_resiude;
	bool is_local_num_in_container = false;
	for (int idx = chain_index_begin; idx < chain_index_end; ++idx) {
		int residue_id = residue_id_container[idx];
		if (residue_id == local_residue_num) {
			is_local_num_in_container = true;
			idx_of_resiude = idx;
			break;
		}
	}
	if (!is_local_num_in_container) {
		std::cerr << "Error: This model does not have the residue." << std::endl;
		std::exit(1);
	}
	
	int result = atom_id_container[idx_of_resiude];
	return result;
}


int cafemol::PDBReader::get_AtomNumof(const std::string& query) {

	// judge whether the string have thein information
	if (chainIDs.empty()) {
		std::cerr << "Error: The string of Chain ID is empty." << std::endl;
		std::exit(1);
	}

	if (query.size() <= 0) {
		std::cerr << "Error: The string of query is empty." << std::endl;
		std::exit(1);
	}

	int result = 0;

	for (int idx = 0; idx < query.size(); ++idx) {
		char chainID = query[idx];

		// search 
		int chain_id_index = chainIDs.find(chainID);
		int chain_index_begin = chain_begin_indices[chain_id_index];
		int chain_index_end = chain_begin_indices[chain_id_index + 1];
		int atom_num_in_chainID = chain_index_end - chain_index_begin;
		result += atom_num_in_chainID;
	
	}
	return result;

}


int cafemol::PDBReader::get_AllAtomNum() {

	int result = get_AtomNumof(chainIDs);
	return result;
}


int cafemol::PDBReader::get_ResidueNumof(const std::string& chainID) {

	if (chainIDs.empty()) {
		std::cerr << "Error: The string of Chain ID is empty." << std::endl;
		std::exit(1);
	}

	int chain_id_index = chainIDs.find(chainID);
	int chain_index_end = chain_begin_indices[chain_id_index + 1] - 1;
	int residue_idx_end = residue_id_container[chain_index_end];
	return residue_idx_end;
}
