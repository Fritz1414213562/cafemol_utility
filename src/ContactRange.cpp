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

	// contact search 
	sout("", "Searching contacted pairs along a DNA chain", "");
	std::vector<std::array<int, 2>> protein_contact_FwRv_ends;
	std::vector<int> dna_contact_atom_Fw_ends;
	std::vector<int> dna_contact_atom_Rv_ends;
	for (const std::array<std::vector<float>, 3>& snapshot : dcd_trajectory) {
		const std::array<std::array<int, 2>, 2>& contact_pairs = contact_search->get_ContactEnds(DNA_IDs, Protein_IDs, snapshot);
		protein_contact_FwRv_ends.push_back({contact_pairs[0][1], contact_pairs[1][1]});
		dna_contact_atom_Fw_ends.push_back(contact_pairs[0][0]);
		dna_contact_atom_Rv_ends.push_back(contact_pairs[1][0]);
	}
	sout(".... Done");

	// convert atom_id of DNA to dsDNA bp of DNA
	sout("Converting the Atom ID to the corresponding Residue ID");
	const std::vector<int>& dna_contact_resi_Fw_ends = topology_reader->convert_DNAID2dsDNAResi(dna_contact_atom_Fw_ends);
	const std::vector<int>& dna_contact_resi_Rv_ends = topology_reader->convert_DNAID2dsDNAResi(dna_contact_atom_Rv_ends);
	sout(".... Done", "", "Calculation ends");
	if ((protein_contact_FwRv_ends.size() != dna_contact_resi_Fw_ends.size()) ||
		(protein_contact_FwRv_ends.size() != dna_contact_resi_Rv_ends.size())) {
		eout("Frame size is not consistent among protein, dna residue ends");
	}

	// output results to file 'output_name'
	const std::size_t& frame_size = protein_contact_FwRv_ends.size();

	sout.output_HyphenBlock("Output the calculation results to " + output_name, BLOCK_SIZE);
	std::ofstream ofs(output_name, std::ios::out);

	ofs << "DNA_ID(fw) PROTEIN_ID(fw) DNA_ID(rv) PROTEIN_ID(rv)" << std::endl;
	
	for (std::size_t iframe = 0; iframe < frame_size; ++iframe) {
		ofs << dna_contact_resi_Fw_ends[iframe] << " " << protein_contact_FwRv_ends[iframe][0] << " "
			<< dna_contact_resi_Rv_ends[iframe] << " " << protein_contact_FwRv_ends[iframe][1] 
			<< std::endl;
	}

	ofs.close();

	sout("Program finished", "");

	return 0;
}
