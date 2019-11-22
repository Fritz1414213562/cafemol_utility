#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/ErrorMessage.hpp"
#include"../include/misc/StandardOutput.hpp"
#include<string>
#include<fstream>
#include<iostream>
#include<memory>
#include<vector>
#include<array>


int main(int argc, char *argv[]) {

    cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
    // command line arguments
    if (argc != 4) error_output("too much or less arguments.");
    cafemol::output_handling::Standard_Output standard_output = cafemol::output_handling::Standard_Output();

    std::array<std::string, 2> input_names = {argv[1], argv[2]};
    std::string output_name = argv[3];

    // default value

    const std::string parameter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/pdns_residue_sort.par";
    const float cutoff = 10.0;
    // block size for standard output
    const std::size_t block_size = 80;

    // suffix
    const std::array<std::string, 2> suffixes = {"dcd", "psf"};
    std::array<std::string, 2> filenames;

    // check each file format
    {
        cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
        filenames = judgement(input_names, suffixes);
    }

    // ---------------------------------------------------------------------------------------------
    // calculation

    // generate instances
    std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
    std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);

    // read protein structure file
    standard_output.output_HyphenBlock("Reading the protein structure information", block_size);

    cafemol::psf_data_type::container_psf_chain_info chain_data = psf_reader->get_ChainInfoOfPSF();
    const std::vector<std::string> atom_name_data = psf_reader->get_AtomNameContainer();
    if (chain_data.empty()) error_output("There is no data in this file.");


    // read the number of DNA chain residues
    standard_output("Reading the beginning and end residue numbers of each DNA chain.");
    std::vector<std::array<int, 2>> dna_chains_start_ends;
    
    int total_DNA_atom_num = 0;
    for (const cafemol::psf_data_type::psf_chain_info& chain_info : chain_data) {
        std::string chain_kind_name = std::get<0>(chain_info);
        if (chain_kind_name == "DNA") {
            std::array<int, 2> chain_start_end = std::get<1>(chain_info);
            dna_chains_start_ends.push_back(chain_start_end);
            total_DNA_atom_num += (chain_start_end[1] - chain_start_end[0] + 1);
        }
    }

    if (dna_chains_start_ends.size() % 2 == 1) error_output("The number of DNA chain is not even.");
    const int dsDNA_number = dna_chains_start_ends.size() / 2;

    std::vector<int> dna_ids;
    for (int dsDNA_idx = 0; dsDNA_idx < dsDNA_number; ++dsDNA_idx) {
        std::array<int, 2> DNA_chainA = dna_chains_start_ends[2 * dsDNA_idx];
        std::array<int, 2> DNA_chainB = dna_chains_start_ends[2 * dsDNA_idx + 1];
        if ((DNA_chainA[1] - DNA_chainA[0]) != (DNA_chainB[1] - DNA_chainB[0])) error_output("This DNA molecules have ssDNA region or The DNA pairing is wrong.");
        for (int dna_id = 0; dna_id < (DNA_chainA[1] - DNA_chainA[0] + 1); ++dna_id) {
            dna_ids.push_back(dna_id + DNA_chainA[0]);
            dna_ids.push_back(DNA_chainB[1] - dna_id);
        }
    }

    standard_output("The number of DNA atoms ->", total_DNA_atom_num);
    standard_output.output_HyphenBlock("", block_size);


    // read the PDNS contact residues
    standard_output.output_HyphenBlock("Reading PDNS contact residue number", block_size);
    std::vector<int> pdns_sorted_residues;

    std::ifstream para_ifs(parameter_name, std::ios::in);
    std::string buffer;
    while (std::getline(para_ifs, buffer)) {
        pdns_sorted_residues.push_back(std::stoi(buffer) + total_DNA_atom_num);
    }
    para_ifs.close();
    if (pdns_sorted_residues.empty()) error_output("No data << ", parameter_name, ">>");

    standard_output("The number of PDNS contacts ->", pdns_sorted_residues.size());
    standard_output.output_HyphenBlock("", block_size);


    // find pdns residues contacting to DNA at each MD steps
    standard_output.output_HyphenBlock("Searching PDNS residues contacting to DNA at each MD steps", block_size);
    std::vector<std::vector<bool>> is_pdns_contacts_vec = dcd_analyzer->is_NativeContacts(pdns_sorted_residues, dna_ids, cutoff);
    standard_output("... Done");
    standard_output.output_HyphenBlock("", block_size);

    standard_output.output_HyphenBlock("Output the results -> " + output_name, block_size);
    standard_output("Mapping PDNS contacts at each MD step");

    std::ofstream ofs(output_name, std::ios::out);

    for (std::size_t iframe = 0; iframe < is_pdns_contacts_vec.size(); ++iframe) {
        for (std::size_t contact_idx = 0; contact_idx < is_pdns_contacts_vec[iframe].size(); ++contact_idx) {
            ofs << iframe << " " << contact_idx + 1 << " " << static_cast<int>(is_pdns_contacts_vec[iframe][contact_idx]) << std::endl;
        }
    }
    ofs.close();

    standard_output("PDNS Mapping is over.");
    standard_output.output_HyphenBlock("", block_size);

    standard_output("Program finished.", "");

    return 0;
}
