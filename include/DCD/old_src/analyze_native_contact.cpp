#include"DCDAnalyzer.hpp"
#include"FileOpenJudge.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<memory>
#include<array>
#include<vector>
#include<numeric>
#include<cmath>
#include<typeinfo>


int main(int argc, char *argv[]) {

    // constant
    const float PI = 4 * atan(1);

    std::string file_name = argv[1];
	std::string output_name = argv[2];
    std::string input_suffix = "dcd";
	// ---------------------------------------------------------------------------------
	// File open
	std::unique_ptr<caf::FileOpenJudge> judgement = std::make_unique<caf::FileOpenJudge>();
	judgement->SuffixJudge(file_name, input_suffix);

	std::string file_prefix = judgement->getFilePrefix();
	// ---------------------------------------------------------------------------------
	// calculation

    std::unique_ptr<caf::DCDAnalyzer> dparser =  std::make_unique<caf::DCDAnalyzer>(file_name);

    std::array<int, 2> dna = {1, 1648}; // The begin and end serial number of dna
    // std::array<int, 2> his8 = {1649, 2628}; // The begin and end serial number of histone octamer
    std::array<int, 2> his8 = {1692, 2628}; // except Histone tail 
	float cut_off = 6.5; // The cutoff of native contact in cafemol

	// remove histone tails
	std::vector<int> erase_nums;
	for (int idx = 1784; idx < 1808; ++idx) {
		erase_nums.push_back(idx);
		erase_nums.push_back(idx + 490);
	}
	for (int idx = 1886; idx < 1901; ++idx) {
		erase_nums.push_back(idx);
		erase_nums.push_back(idx + 490);
	}
	for (int idx = 2001; idx < 2050; ++idx) {
		erase_nums.push_back(idx);
		erase_nums.push_back(idx + 490);
	}
	for (int idx = 2139; idx < 2183; ++idx) erase_nums.push_back(idx);


	std::vector<int> nat_contacts = dparser->searchNativeContact(dna, his8, erase_nums, cut_off);
	int frame_size = nat_contacts.size();
	std::ofstream ofs(output_name, std::ios::out);

	for (std::size_t frame_i = 0; frame_i < frame_size; ++frame_i) {
		ofs << nat_contacts[frame_i] << std::endl;
	}

    return 0;
}

