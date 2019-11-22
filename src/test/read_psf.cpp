#include"../../include/PSF/PSFReader.hpp"
#include<string>
#include<iostream>
#include<memory>
#include<vector>
#include<string>
#include<array>


int main(int argc, char* argv[]) {

	std::string psf = argv[1];

	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(psf);

//	psf_reader->read_ATOM_BOND();
//	std::vector<std::string> res = psf_reader->search_ChainKind();

	std::vector<cafemol::psf_chain_info> result = psf_reader->get_ChainInfoOfPSF();

	for (const cafemol::psf_chain_info& pci : result) {
		std::string chain_name = std::get<0>(pci);
		std::array<int, 2> chain_start_end = std::get<1>(pci);
		std::cout << chain_name << " " << chain_start_end[0] << " " << chain_start_end[1] << std::endl;
	}
	
	return 0;
}
