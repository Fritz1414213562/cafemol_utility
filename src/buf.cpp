#include"../include/NINFO/NinfoReader.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/PDB/PDBReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/CafemolLibrary.hpp"
#include<array>
#include<vector>
#include<fstream>
#include<memory>
#include<iostream>
#include<set>


int main(int argc, char *argv[]) {

	// for file stream
	std::string ninfo_name = argv[1];
	std::string psf_name = argv[2];
	std::string pdb_name = argv[3];
	std::string output_name = argv[4];

	// contact
	std::string ninfo_suffix = "ninfo";
	std::string psf_suffix = "psf";
	std::string pdb_suffix = "pdb";

	// -----------------------------------------------------------------------------
	// File open
	{
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(ninfo_name, ninfo_suffix);
		judgement->SuffixJudge(psf_name, psf_suffix);
		judgement->SuffixJudge(pdb_name, pdb_suffix);
	}

	// -----------------------------------------------------------------------------
	// calculation

	// generate instance
	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(ninfo_name);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(psf_name);
	std::unique_ptr<cafemol::PDBReader> pdb_reader = std::make_unique<cafemol::PDBReader>(pdb_name);

	// read file data
	cafemol::ninfo_data_type::container_pdns pdns_data = ninfo_reader->read_NinfoPdns();
	std::vector<cafemol::psf_data_type::psf_chain_info> all_chains = psf_reader->get_ChainInfoOfPSF();
	pdb_reader->read_PDBData();


	std::vector<std::vector<int>> pdns_histone_chains_ids;
	std::vector<int> pdns_histone_ids;
	int current_chain = 0;
	for (const cafemol::ninfo_data_type::pdns_tuple& pdns_line_data : pdns_data) {
		std::array<int, 2> pdns_mp = std::get<2>(pdns_line_data.line_data);
		int pdns_chain_num = std::get<1>(pdns_line_data.line_data)[0]
		if (current_chain == 0) {
			current_chain = pdns_chain_num;
		}
		else if (pdns_chain_num != current_chain) {
			pdns_histone_chain_ids.push_back(pdns_histone_ids);
			current_chain = pdns_chain_num;
			pdns_histone_ids.clear();
		}
		else continue;
	}

	return 0;
}
