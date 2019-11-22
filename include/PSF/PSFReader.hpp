#ifndef PSF_READER_HPP
#define PSF_READER_HPP
#include"PSFParser.hpp"
#include"../misc/ErrorMessage.hpp"
#include<string>
#include<vector>
#include<algorithm>
#include<iostream>
#include<tuple>


namespace cafemol {

namespace psf_data_type {
	using chain_info_tuple = std::tuple<std::string, std::array<int, 2>>;
	using psf_chain_info = chain_info_tuple;
	using container_psf_chain_info = std::vector<psf_chain_info>;

}

class PSFReader : public PSFParser {

public:

	PSFReader() = default;
	PSFReader(const std::string& filename);
	
	std::vector<psf_data_type::psf_chain_info> get_ChainInfoOfPSF();
	std::vector<std::array<int, 2>> get_ResidueBeginEnds();
	int get_ResidueNumof(const int& query_atom_id);
	int get_LocalResidueNumof(const int& query_atom_id);
	std::vector<std::string> get_AtomNameContainer();
	std::vector<int> convert_DNAID2dsDNAResi(const std::vector<int>& dna_ids);

// semi-public
	void read_ATOM_BOND();
	std::vector<std::string> search_ChainKind();

//// If you want to use this method, you must uncomment in this and cpp file.
//	void test_ChainSearch();

private:
	// !NATOM
	std::vector<int> atom_id_container;
	std::vector<std::string> segment_name_container;
	std::vector<int> residue_num_container;
	std::vector<std::string> residue_name_container;
	std::vector<std::string> atom_name_container;
	std::vector<std::string> atom_type_container;
	std::vector<std::string> charge_container;
	std::vector<std::string> mass_container;

	// !NBOND
	std::vector<std::array<int, 2>> bond_pairs;

	// about molecule
	// store the INDEX of "atom_id_container", not the ATOM ID.
	std::vector<std::array<int, 2>> dna_chains_start_ends;
	std::vector<std::array<int, 2>> pro_chains_start_ends;
	// store the ATOM ID.
	std::vector<std::array<int, 2>> chain_start_end_ids;
    // store the ATOM Name
    std::vector<std::string> chain_names;

	void read_AtomInfo(const int& atom_num);
	void read_BondInfo(const int& bond_num);

	void get_BondPair(const std::string& line);
	void skip_RowsUntil(const std::string& block_name, int& block_num);

	void search_DNAAtoms();
	void search_ProteinAtoms();

	void search_Chain();
	void judge_ChainKind(std::vector<std::string>& chain_name_container, const std::vector<std::array<int, 2>>& target_chains_start_ends, const std::string& chain_name);

	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
};
}

#endif /* PSF_READER_HPP */
