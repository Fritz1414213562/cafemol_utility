#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PDB/PDBReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<memory>
#include<array>
#include<vector>


int main(int argc, char *argv[]) {

	if (argc != 5) {
		std::cerr << "too much or less arguments" << std::endl;
		std::exit(1);
	}

	// for file stream
    std::string dcd_name = argv[1];
	std::string pdb_name = argv[2];
	std::string output_name = argv[3];
	float cutoff = std::stof(argv[4]);
    std::string dcd_suffix = "dcd";
	std::string pdb_suffix = "pdb";
	// ---------------------------------------------------------------------------------
	// local constant
//	float cutoff = 6.5; // The cutoff length of native contact in cafemol
//	float cutoff = 8.0; // The cutoff length of native contact in cafemol
	// local ID in a Histone monomer.
	int local_His3_begin = 45;
	int local_His3_end = 135;
	int local_His4_begin = 31;
	int local_His4_end = 102;
	int local_His2A_begin = 16;
	int local_His2A_end = 103;
	int local_His2B_begin = 37;
	int local_His2B_end = 125;

	// chain ID
	std::string His8ChainIDs = "CDEFGHIJ";
	std::string His4merIDs = "CDEF";
	std::string His3ID = "C";
	std::string His4ID = "D";
	std::string His2AID = "E";
	std::string His2BID = "F";
	std::string DNAChainID = "AB";
	std::string DNAChainIDA = "A";

	// ---------------------------------------------------------------------------------
	// File open
	{
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(dcd_name, dcd_suffix);
		judgement->SuffixJudge(pdb_name, pdb_suffix);
	}

	// ---------------------------------------------------------------------------------
	// calculation

	// generate instance
    std::unique_ptr<cafemol::DCDAnalyzer> dparser =  std::make_unique<cafemol::DCDAnalyzer>(dcd_name);
	std::unique_ptr<cafemol::PDBReader> pparser = std::make_unique<cafemol::PDBReader>(pdb_name);
	pparser->read_PDBData();

	std::cout << "-------Reading the information on DNA from PDB------------------------------------" << std::endl;
	int dna_atom_size = pparser->get_AtomNumof(DNAChainID);
	std::cout << "The number of atoms belonging to dsDNA -> " << dna_atom_size << std::endl;
	int dna_A_size = pparser->get_AtomNumof(DNAChainIDA);
	int dna_resi_size = pparser->get_ResidueNumof(DNAChainIDA);
	std::cout << "The length of base pair -> " << dna_resi_size << std::endl;

	std::vector<int> dna_ids(dna_atom_size);
	dna_ids[0] = 1;
	dna_ids[1] = 2;
	dna_ids[dna_atom_size - 2] = dna_A_size + 2;
	dna_ids[dna_atom_size - 1] = dna_A_size + 1;

	for (int idx = 0; idx < (dna_atom_size - 4) / 6; ++idx) {
		dna_ids[6 * idx + 2] = dna_atom_size - 3 * idx; 
		dna_ids[6 * idx + 3] = dna_atom_size - 3 * idx - 1;
		dna_ids[6 * idx + 4] = dna_atom_size - 3 * idx - 2;
		dna_ids[6 * idx + 5] = 3 * idx + 3;
		dna_ids[6 * idx + 6] = 3 * idx + 4;
		dna_ids[6 * idx + 7] = 3 * idx + 5;
	}

	std::cout << "-------Reading the information on Histone 8mer------------------------------------" << std::endl;
	int His8_atom_begin = dna_atom_size + 1;
	int His8_atom_end = pparser->get_AllAtomNum();
	std::cout << "The range of Histone octamer ID -> " << His8_atom_begin << "-" << His8_atom_end << std::endl;

	// store histone ids removing tails' region
	int His3_begin = pparser->get_AtomIDof(local_His3_begin, His3ID);
	int His3_end = pparser->get_AtomIDof(local_His3_end, His3ID);
	int His4_begin = pparser->get_AtomIDof(local_His4_begin, His4ID);
	int His4_end = pparser->get_AtomIDof(local_His4_end, His4ID);
	int His2A_begin = pparser->get_AtomIDof(local_His2A_begin, His2AID);
	int His2A_end = pparser->get_AtomIDof(local_His2A_end, His2AID);
	int His2B_begin = pparser->get_AtomIDof(local_His2B_begin, His2BID);
	int His2B_end = pparser->get_AtomIDof(local_His2B_end, His2BID);

	// histone 4mer's number
	int His4mer_num = pparser->get_AtomNumof(His4merIDs);

	std::vector<int> His8_ids;
//	for (int id = His8_atom_begin; id <= His8_atom_end; ++id) {
//		His8_ids.push_back(id);
//	}

	std::cout << "The range of Histone ID 'C' -> " << His3_begin << "-" << His3_end << std::endl;
	std::cout << "The range of Histone ID 'D' -> " << His4_begin << "-" << His4_end << std::endl;
	std::cout << "The range of Histone ID 'E' -> " << His2A_begin << "-" << His2A_end << std::endl;
	std::cout << "The range of Histone ID 'F' -> " << His2B_begin << "-" << His2B_end << std::endl;
	std::cout << "The range of Histone ID 'G' -> " << His3_begin + His4mer_num << "-" << His3_end + His4mer_num << std::endl;
	std::cout << "The range of Histone ID 'H' -> " << His4_begin + His4mer_num << "-" << His4_end + His4mer_num << std::endl;
	std::cout << "The range of Histone ID 'I' -> " << His2A_begin + His4mer_num << "-" << His2A_end + His4mer_num << std::endl;
	std::cout << "The range of Histone ID 'J' -> " << His2B_begin + His4mer_num << "-" << His2B_end + His4mer_num << std::endl;

	for (int id = His3_begin; id <= His3_end; ++id) {
		His8_ids.push_back(id);
		His8_ids.push_back(id + His4mer_num);
	}
	
	for (int id = His4_begin; id <= His4_end; ++id) {
		His8_ids.push_back(id);
		His8_ids.push_back(id + His4mer_num);
	}
	
	for (int id = His2A_begin; id <= His2A_end; ++id) {
		His8_ids.push_back(id);
		His8_ids.push_back(id + His4mer_num);
	}
	
	for (int id = His2B_begin; id <= His2B_end; ++id) {
		His8_ids.push_back(id);
		His8_ids.push_back(id + His4mer_num);
	}
	
	std::cout << "-------Calculating the 1st base pair contacted to Histone octamer-----------------" << std::endl;
	std::vector<std::array<int, 2>> closest_pairs = dparser->get_1stNativeContactedResidue(dna_ids, His8_ids, cutoff);
	std::vector<int> closest_ids;
	for (const std::array<int, 2>& closest_pair : closest_pairs) {
		int closest_dna_id = closest_pair[0];
		closest_ids.push_back(closest_dna_id);
	}
	//std::cout << closest_ids.size() << std::endl;
	std::cout << "-------Caluculation is over-------------------------------------------------------" << std::endl;
	std::cout << "-------Checking the number of atoms in each file. --------------------------------" << std::endl;
	// check atom number written in each file
	int dcd_atom = dparser->get_AtomNumber();
	int pdb_atom = pparser->get_AllAtomNum();
	std::cout << "The number of atoms in PDB file -> " << pdb_atom << std::endl;
	std::cout << "The number of atoms in DCD file -> " << dcd_atom << std::endl;

	if (dcd_atom != pdb_atom) {
		std::cerr << "Error: The number of atoms in PDB is not consistent with that in DCD." << std::endl;
		std::exit(1);
	}

	std::cout << "-------Making an output file------------------------------------------------------" << std::endl;

	std::ofstream ofs(output_name, std::ios::out);

	for (const int& closest_id : closest_ids) {
		int closest_bp;
		if (closest_id > dna_A_size) {
			closest_bp = pparser->get_ResidueIDof(closest_id);
			closest_bp = dna_resi_size - closest_bp + 1;
			ofs << closest_bp << std::endl;
		}
		else if ((closest_id <= dna_A_size) && (0 < closest_id)) {
			closest_bp = pparser->get_ResidueIDof(closest_id);
			ofs << closest_bp << std::endl;
		}

		else if (closest_id <= 0) {
			ofs << "nan" << std::endl;
		}
	}

	ofs.close();
	std::cout << "The program finished" << std::endl;
	std::cout << std::endl;

    return 0;
}

