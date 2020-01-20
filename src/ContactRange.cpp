#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<DCD/DCDReader.hpp>
#include<PSF/PSFReader.hpp>
#include<detectContacts.hpp>
#include<FileOpenJudge.hpp>
#include<CafemolName.hpp>
#include<fstream>
#include<string>
#include<memory>
#include<vector>
#include<array>


int main(int argc, char *argv[]) {

	// ------------------------------------------------------------------------------------------
	// command line arguments

	const int& max_argc = 5;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::array<std::string, 2>& input_names = {argv[1], argv[2]};
	const std::string& output_name = argv[3];
	const float& cutoff_len = std::stof(argv[4]);

	// local variables
	const std::size_t& BLOCK_SIZE = 90;

	const std::array<std::string, 2>& suffixes = {"dcd", "psf"};
	std::array<std::string, 2> filenames;

	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(input_names, suffixes);
	}

	// ------------------------------------------------------------------------------------------
	// calculation
	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Calculation starts.", "");

	// trajectory
	std::unique_ptr<cafemol::DCDReader> trajectory_reader = std::make_unique<cafemol::DCDReader>(filenames[0]);
	// topology
	std::unique_ptr<cafemol::PSFReader> topology_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);
	// contact search
	std::unique_ptr<cafemol::ContactSearcher> contact_search = std::make_unique<cafemol::ContactSearcher>(cutoff_len);

	// reading a trajectory
	sout("Reading a trajectory", "");
	const cafemol::Trajectory& dcd_trajectory = trajectory_reader->get_Trajectory();
	sout(".... Done", "");

	// reading a topology
	sout("Reading the topology", "");
	const std::vector<std::size_t>& DNA_IDs = topology_reader->get_AlignedDNAIDs();
	const std::vector<std::size_t>& Protein_IDs = topology_reader->get_ProteinIDs();
	sout(".... Done", "");

	// contact search and output data
	sout.output_HyphenBlock("Output the calculation results to " + output_name, BLOCK_SIZE);
	sout("", "Searching contacted pairs along a DNA chain", "");
	std::ofstream ofs(output_name, std::ios::out);

	ofs << "DNA_ID(fw) PROTEIN_ID(fw) DNA_ID(rv) PROTEIN_ID(rv)" << std::endl;
	for (const std::array<std::vector<float>, 3>& snapshot : dcd_trajectory) {
		const std::array<std::array<int, 2>, 2>& contact_pairs = contact_search->get_ContactEnds(DNA_IDs, Protein_IDs, snapshot);
		ofs << contact_pairs[0][0] << " " << contact_pairs[0][1] << " " << contact_pairs[1][0] << " " << contact_pairs[1][1] << std::endl;
	}
	ofs.close();

	sout("Program finished", "");

	return 0;
}
