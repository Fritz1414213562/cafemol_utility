#include"PSFReader.hpp"


// -------------------------------------------------------------------------------------------
// public


cafemol::PSFReader::PSFReader(const std::string& filename) : PSFParser(filename) {}

void cafemol::PSFReader::read_ATOM_BOND() {

	if (!ifs.is_open()) {
		std::cerr << "Error: The file is not opened." << std::endl;
		std::exit(1);
	}

	int atom_num;
	skip_RowsUntil("!NATOM", atom_num);
	read_AtomInfo(atom_num);

	int bond_num;
	skip_RowsUntil("!NBOND", bond_num);
	read_BondInfo(bond_num);
}


std::vector<std::string> cafemol::PSFReader::search_ChainKind() {

	// search chain kind
	search_DNAAtoms();
	search_ProteinAtoms();
	search_Chain();

	//for (const std::array<int, 2>& arr : dna_chains_start_ends) {
	//	std::cout << arr[0] << " " << arr[1] << std::endl;
	//}
	//std::cout <<  pro_chains_start_ends.size() << std::endl;
	//for (const std::array<int, 2>& arr : pro_chains_start_ends) {
	//	if (arr.empty()) {
	//		std::cout << "empty" << std::endl;
	//		break;
	//	}
	//	std::cout << arr[0] << " " << arr[1] << std::endl;
	//}

	std::vector<std::string> result(chain_start_end_ids.size());
	judge_ChainKind(result, dna_chains_start_ends, "DNA");
	judge_ChainKind(result, pro_chains_start_ends, "Protein");
	//for (const std::string& Str : result) {
	//	if (Str.empty()) {
	//		std::cout << "Empty" << std::endl;
	//	}
	//	else {
	//		std::cout << Str << std::endl;
	//	}
	//}

	return result;
}


std::vector<cafemol::psf_data_type::psf_chain_info> cafemol::PSFReader::get_ChainInfoOfPSF() {

	read_ATOM_BOND();
	chain_names = search_ChainKind();

	if (chain_start_end_ids.size() != chain_names.size()) {
		std::cerr << "Error: Something wrong! The size of Chain IDs is not consistent with that of Chain names" << std::endl;
		std::exit(1);
	}

	std::vector<cafemol::psf_data_type::psf_chain_info> result(chain_start_end_ids.size());

	for (std::size_t idx = 0; idx < chain_start_end_ids.size(); ++idx) {
		std::get<0>(result[idx]) = chain_names[idx];
		std::get<1>(result[idx]) = chain_start_end_ids[idx];
	}

	return result;
}


int cafemol::PSFReader::get_ResidueNumof(const int& query_atom_id) {

	if (atom_id_container.empty() || residue_num_container.empty()) {
		std::cerr << "Error: have not yet read the file." << std::endl;
		std::exit(1);
	}
	
	int idx = 0;
	for (const int& atom_id : atom_id_container) {
		if (atom_id == query_atom_id) break;
		++idx;
	}
	if (idx >= atom_id_container.size()) {
		std::cerr << "Error: The query ID '" << query_atom_id << "' could not be found." << std::endl;
		std::exit(1);
	}
	int result = residue_num_container[idx];
	return result;
}


int cafemol::PSFReader::get_LocalResidueNumof(const int& query_atom_id) {

	int global_residue_num = get_ResidueNumof(query_atom_id);

	int residue_idx = 0;
	for (const std::array<int, 2>& chain_start_end : chain_start_end_ids) {
		if ((query_atom_id >= chain_start_end[0]) && (query_atom_id <= chain_start_end[1])) break;
		++residue_idx;
	}
	int residue_begin_id = get_ResidueNumof(chain_start_end_ids[residue_idx][0]);

	int local_residue_num = global_residue_num - residue_begin_id + 1;
	if (local_residue_num < 0) {
		std::cerr << "Error: The residue search did not work completely." << std::endl;
		std::exit(1);
	}

	return local_residue_num;
}


std::vector<std::array<int, 2>> cafemol::PSFReader::get_ResidueBeginEnds() {

	if (chain_start_end_ids.empty()) {
		std::cerr << "Error: have not yet read the file." << std::endl;
		std::exit(1);
	}
	std::vector<std::array<int, 2>> result;

	for (const std::array<int, 2>& chain_start_end : chain_start_end_ids) {
		int residue_start = get_ResidueNumof(chain_start_end[0]);
		int residue_end = get_ResidueNumof(chain_start_end[1]);
		result.push_back({residue_start, residue_end});
	}
	return result;
}


std::vector<std::string> cafemol::PSFReader::get_AtomNameContainer() {
	if (chain_start_end_ids.empty()) error_output("have not yet read the file.");
	return atom_name_container;
}


std::vector<int> cafemol::PSFReader::convert_DNAID2dsDNAResi(const std::vector<int>& dna_ids) {
	if (chain_start_end_ids.empty()) error_output("have not yet read the file.");
    else if (chain_names.empty()) error_output("have not yet read chain names in the file.");

    std::vector<std::array<int, 2>> residue_num_begin_end = get_ResidueBeginEnds();
	std::vector<std::array<int, 2>> dna_start_ends;
    const std::size_t chain_size = chain_names.size();
	for (std::size_t chain_idx = 0; chain_idx < chain_size; ++chain_idx) {
		if (chain_names[chain_idx] == "DNA") {
			std::array<int, 2> chain_start_end = chain_start_end_ids[chain_idx];
			dna_start_ends.push_back(chain_start_end);
		}
	}
	if (dna_start_ends.empty()) error_output("This model have no DNA molecule.");
	
	std::vector<int> result;
	for (const int& dna_id : dna_ids) {
		if (dna_id <= 0) {
			result.push_back(-1);
			continue;
		}

		int dna_bp;
		int DNA_chain_num = 0;
		int dsDNA_length = 0;
		bool is_DNA = false;
		
		// determine which chain DNA atom id belongs to.
		for (const std::array<int, 2>& dna_start_end : dna_start_ends) {
            if ((dna_id >= dna_start_end[0]) && (dna_id <= dna_start_end[1])) {
                dsDNA_length = residue_num_begin_end[DNA_chain_num][1] - residue_num_begin_end[DNA_chain_num][0] + 1;
                ++DNA_chain_num;
                is_DNA = true;
                break;
            }
            ++DNA_chain_num;
        }
        if (!is_DNA) {
            std::cerr << "Error: not found the output DNA residue number." << std::endl;
            std::exit(1);
        }

        if (DNA_chain_num % 2 == 0) {
            dna_bp = get_LocalResidueNumof(dna_id);
            dna_bp = dsDNA_length - dna_bp + 1;
        }
        else if (DNA_chain_num % 2 == 1) {
            dna_bp = get_LocalResidueNumof(dna_id);
        }
        else {
            std::cerr << "Error: Something wrong." << std::endl;
            std::exit(1);
        }
        result.push_back(dna_bp);
    }
    return result;
}


//void cafemol::PSFReader::test_ChainSearch() {
//	read_ATOM_BOND();
//	search_ChainKind();
//	search_Chain();
//
//	for (const std::array<int, 2>& chain_start_end : chain_start_end_ids) {
//		std::cout << "write chain ids of its begin and end." << std::endl;
//		std::cout << chain_start_end[0] << " " << chain_start_end[1] << std::endl;
//	}
//}
//
// -------------------------------------------------------------------------------------------
// private


void cafemol::PSFReader::read_AtomInfo(const int& atom_num) {

	atom_id_container.resize(atom_num);
	segment_name_container.resize(atom_num);
	residue_num_container.resize(atom_num);
	residue_name_container.resize(atom_num);
	atom_name_container.resize(atom_num);
	atom_type_container.resize(atom_num);
	charge_container.resize(atom_num);
	mass_container.resize(atom_num);

	for (std::size_t line_idx = 0; line_idx < atom_num; ++line_idx) {
		std::string line;
		std::getline(ifs, line);
		atom_id_container[line_idx] = get_AtomID(line);
		segment_name_container[line_idx] = get_SegmentName(line);
		residue_num_container[line_idx] = get_ResidueNum(line);
		residue_name_container[line_idx] = get_ResidueName(line);
		atom_name_container[line_idx] = get_AtomName(line);
		atom_type_container[line_idx] = get_AtomType(line);
		charge_container[line_idx] = get_Charge(line);
		mass_container[line_idx] = get_Mass(line);
	}
}


void cafemol::PSFReader::read_BondInfo(const int& bond_num) {
	int row_size = (bond_num / 4) + 1;

	for (std::size_t line_idx = 0; line_idx < row_size; ++line_idx) {
		std::string line;
		std::getline(ifs, line);
		get_BondPair(line);
	}
}


void cafemol::PSFReader::get_BondPair(const std::string& line) {
	int line_length = line.size();
	line_length = (line_length / 8) * 8;

	for (std::size_t idx = 0; idx < line_length; idx += 16) {
		std::string buf1 = line.substr(idx, 8);
		std::string buf2 = line.substr(idx + 8, 8);
		bond_pairs.push_back({std::stoi(buf1), std::stoi(buf2)});
	}
}


void cafemol::PSFReader::skip_RowsUntil(const std::string& block_name, int& block_num) {

	std::string line;
	while (std::getline(ifs, line)) {

		if (line.empty()) continue;
		std::string block_kind = get_BlockKind(line);

		if (block_kind[0] != '!') continue;
		else if (block_kind != block_name) {
			std::cerr << "Error: This file does not have '" << block_name << "' row" << std::endl;
			std::exit(1);
		}
		else if (block_kind == block_name) {
			block_num = get_BlockNum(line);
			break;
		}
		else {
			std::cerr << "Error: Something wrong. Maybe, block_name '" << block_name << "' is wrong." << std::endl;
			std::exit(1);
		}
	}
}


void cafemol::PSFReader::search_DNAAtoms() {

	if (atom_name_container.empty()) {
		std::cerr << "Error: You don't read the PSF file or ATOM info." << std::endl;
		std::exit(1);
	}

	dna_chains_start_ends.clear();

	std::size_t atom_num = atom_name_container.size();
	bool is_on_DNA = false;
	std::array<int, 2> dna_start_end = {0, 0};
	int previous_atom_idx = 0;

	for (std::size_t idx = 0; idx < atom_num; ++idx) {
		std::string atom_name = atom_name_container[idx];
		int atom_idx = idx;
		if (((atom_name == "DP") 
		 || (atom_name == "DS") 
		 || (atom_name == "DB"))
		 && (is_on_DNA) && (idx < atom_num - 1)) {
			previous_atom_idx = atom_idx;
			continue;
		}
		else if (((atom_name == "DP") 
			   || (atom_name == "DS") 
			   || (atom_name == "DB"))
			   && (is_on_DNA) && (idx == atom_num - 1)) {
			dna_start_end[1] = atom_idx;
			dna_chains_start_ends.push_back(dna_start_end);
		}
		else if (((atom_name == "DP")
			  || (atom_name == "DS")
			  || (atom_name == "DB"))
			  && (!is_on_DNA)) {
			is_on_DNA = true;
			dna_start_end[0] = atom_idx;
			continue;
		}
		else if (is_on_DNA) {
			is_on_DNA = false;
			dna_start_end[1] = previous_atom_idx;
			dna_chains_start_ends.push_back(dna_start_end);
			dna_start_end = {0, 0};
			continue;
		}
		else continue;
	}
}


void cafemol::PSFReader::search_ProteinAtoms() {

	if (atom_name_container.empty()) {
		std::cerr << "Error: You don't read the PSF file or ATOM info." << std::endl;
		std::exit(1);
	}
	
	pro_chains_start_ends.clear();

	std::size_t atom_num = atom_name_container.size();
	bool is_on_Protein = false;
	std::array<int, 2> pro_start_end = {0, 0};
	int previous_atom_idx = 0;

	for (std::size_t idx = 0; idx < atom_num; ++idx) {
		std::string atom_name = atom_name_container[idx];
		int atom_idx = idx;
		if ((atom_name == "CA") && (is_on_Protein) && (idx < atom_num - 1)) {
			previous_atom_idx = atom_idx;
			continue;
		}
		else if ((atom_name == "CA") && (is_on_Protein) && (idx == atom_num - 1)) {
			pro_start_end[1] = atom_idx;
			pro_chains_start_ends.push_back(pro_start_end);
		}
		else if ((atom_name == "CA") && (!is_on_Protein)) {
			is_on_Protein = true;
			pro_start_end[0] = atom_idx;
			continue;
		}
		else if (is_on_Protein) {
			is_on_Protein = false;
			pro_start_end[1] = previous_atom_idx;
			pro_chains_start_ends.push_back(pro_start_end);
			pro_start_end = {0, 0};
			continue;
		}
		else continue;
	}
}


void cafemol::PSFReader::search_Chain() {

	if (bond_pairs.empty()) {
		std::cerr << "Error: You don't read !NBOND block." << std::endl;
		std::exit(1);
	}

	if (!chain_start_end_ids.empty()) {
		chain_start_end_ids.clear();
	}

	std::array<int, 2> chain_start_end = {std::min(bond_pairs[0][0], bond_pairs[0][1]), 0};
	std::array<int, 2> previous_bond_pair = bond_pairs[0];
	std::size_t bond_pairs_num = bond_pairs.size();
	for (std::size_t idx = 1; idx < bond_pairs_num; ++idx) {
		if ((bond_pairs[idx][0] == previous_bond_pair[0]) ||
			(bond_pairs[idx][0] == previous_bond_pair[1]) ||
			(bond_pairs[idx][1] == previous_bond_pair[0]) ||
			(bond_pairs[idx][1] == previous_bond_pair[1])) {
			previous_bond_pair = bond_pairs[idx];
			continue;
		}
		else {
			chain_start_end[1] = std::max(previous_bond_pair[0], previous_bond_pair[1]);
			chain_start_end_ids.push_back(chain_start_end);
			chain_start_end[0] = std::min(bond_pairs[idx][0], bond_pairs[idx][1]);
			previous_bond_pair = bond_pairs[idx];
			continue;
		}
	}

	chain_start_end[1] = std::max(previous_bond_pair[0], previous_bond_pair[1]);
	chain_start_end_ids.push_back(chain_start_end);
}


void cafemol::PSFReader::judge_ChainKind(std::vector<std::string>& chain_name_container, const std::vector<std::array<int, 2>>& target_chains_start_ends, const std::string& chain_name) {


	for (std::size_t groups_count = 0; groups_count < target_chains_start_ends.size(); ++groups_count) {
		std::size_t chain_idx = 0;
		for (const std::array<int, 2>& chain_start_end : chain_start_end_ids) {
			std::array<int, 2> target_chain_start_end = target_chains_start_ends[groups_count];
			int target_chain_start = atom_id_container[target_chain_start_end[0]];
			int target_chain_end = atom_id_container[target_chain_start_end[1]];
			// std::cout << chain_idx << std::endl;
	
			if ((chain_start_end[0] >= target_chain_start) && (chain_start_end[1] <= target_chain_end) && (chain_name_container[chain_idx].empty())) {
				chain_name_container[chain_idx] = chain_name;
				++chain_idx;
			}

			else if ((chain_start_end[0] > target_chain_end) || (chain_start_end[1] < target_chain_start)) {
				++chain_idx;
				continue;
			}

			else if (!chain_name_container[chain_idx].empty()) {
				std::cerr << "Error: The chain group is collide with the previous one." << std::endl;
				std::exit(1);
			}

			else {
				std::cerr << "Error: DNA groups lie between two different chains." << std::endl;
				std::exit(1);
			}
	
		}
	}

}


