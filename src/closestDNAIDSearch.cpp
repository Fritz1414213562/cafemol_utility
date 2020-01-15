#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<DCD/DCDReader.hpp>
#include<PSF/PSFReader.hpp>
#include<detectContacts.hpp>
#include<FileOpenJudge.hpp>
#include<CafemolName.hpp>
#include<array>
#include<vector>
#include<string>
#include<memory>
#include<fstream>


int main(int argc, char *argv[]) {

	const int& max_argc = 6;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::array<std::string, 2>& inputnames = {argv[1], argv[2]};
	const std::string& output_name = argv[3];
	const std::size_t& atom_id = std::stoi(argv[4]);
	const float& cutoff_len = std::stof(argv[5]);

	const std::size_t& BLOCK_SIZE = 90;

	const std::array<std::string, 2>& suffixes = {"dcd", "psf"};
	std::array<std::string, 2> filenames;
	
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(inputnames, suffixes);
	}

	// -----------------------------------------------------------------------------------------
	// calculation
	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Calculation starts.", "");

	// trajectory
	std::unique_ptr<cafemol::DCDReader> trajectory_reader = std::make_unique<cafemol::DCDReader>(filenames[0]);
	// topology
	std::unique_ptr<cafemol::PSFReader> topology_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);
	// contact search
	std::unique_ptr<cafemol::ContactSearcher> contact_search = std::make_unique<cafemol::ContactSearcher>(cutoff_len);


	// read trajectory
	sout("Reading a trajectory");
	const cafemol::Trajectory& dcd_trajectory = trajectory_reader->get_Trajectory();
	sout(".... Done", "");

	// read topology and DNA IDs
	sout("Reading a Topology");
	const std::vector<std::size_t>& DNA_IDs = topology_reader->get_DNAIDs();
	sout(".... Done", "");

	// contact search
	sout("Searching the contacts closest to AtomID:" + std::to_string(atom_id) + ".");
	const std::vector<int>& closest_atomid_list = contact_search->get_TimeSeriesOfClosestID2(atom_id, DNA_IDs, dcd_trajectory);
	sout(".... Done", "");

	// convert Atom IDs to Residue ID
	sout("Converting the Atom ID to the corresponding Residue ID");
	const std::vector<int>& closest_resi_list = topology_reader->convert_DNAID2dsDNAResi(closest_atomid_list);	
	sout(".... Done", "", "Calculation ends");

	// output the result to file
	sout.output_HyphenBlock("Output the calculation results to " + output_name, BLOCK_SIZE);
	std::ofstream ofs(output_name, std::ios::out);

	for (const int& closest_resi : closest_resi_list) {
		if (closest_resi <= 0) ofs << "nan" << std::endl;
		else ofs << closest_resi << std::endl;
	}
//	for (const int& closest_atomid : closest_atomid_list) {
//		if (closest_atomid <= 0) ofs << "nan" << std::endl;
//		else ofs << closest_atomid << std::endl;
//	}
	ofs.close();
	sout("", ".... Done", "");
	sout("Program finished", "");

	return 0;
}
