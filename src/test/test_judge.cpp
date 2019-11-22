#include"../../include/misc/FileOpenJudge.hpp"
#include<iostream>
#include<array>
#include<string>
#include<memory>

int main(int argc, char* argv[]) {

	cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
//	std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
//	std::array<std::string, 4> filenames = {"test.dcd",
//											"test.psf",
//											"test.pdb",
//											"test.tes"};
//
	std::array<std::string, 4> filenames = {argv[1],
											argv[2],
											argv[3],
											argv[4]};
	std::array<std::string, 4> suffixes = {"dcd",
										   "pdb",
										   "tes",
										   "psf"};
	
//	std::array<std::string, 4> res = judgement.SuffixesJudge(filenames, suffixes);
	std::array<std::string, 4> res = judgement(filenames, suffixes);
//	judgement->SuffixesJudge(filenames, suffixes);
	for (const std::string& res_str : res) {
		std::cout << res_str << std::endl;
	}

	return 0;
}
