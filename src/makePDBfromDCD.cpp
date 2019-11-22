#include"../include/DCD/DCDReader.hpp"
#include"../include/PDB/PDBWriter.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include<string>
#include<iostream>
#include<memory>
#include<array>
#include<vector>


int main(int argc, char *argv[]) {

	// command line argument
	std::string dcd_name = argv[1];
	std::string pdb_name = argv[2];
	std::string output_name = argv[3];
	std::string s_frame = argv[4];
	int i_frame = stoi(s_frame);

	// -----------------------------------------------------------------------------------
	// for check file format
	{
		std::string dcd_suffix = "dcd";
		std::string pdb_suffix = "pdb";
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(dcd_name, dcd_suffix);
		judgement->SuffixJudge(pdb_name, pdb_suffix);
	}

	// -----------------------------------------------------------------------------------
	// make PDB file
	
	// generate instances
	std::unique_ptr<cafemol::DCDReader> dcdreader = std::make_unique<cafemol::DCDReader>(dcd_name);
	std::unique_ptr<cafemol::PDBWriter> pdbwriter = std::make_unique<cafemol::PDBWriter>();

	// program start
	std::cout << "----------------------- Taking a snapshot from DCD file --------------------------" << std::endl;
	std::array<std::vector<float>, 3> snapshot = dcdreader->get_SnapShot_at(i_frame);
	std::cout << "----------------------- Making a PDB file from template --------------------------" << std::endl;
	pdbwriter->makePDBFileFrom(pdb_name, output_name, snapshot);
	std::cout << "The program finished" << std::endl;

	return 0;
}
