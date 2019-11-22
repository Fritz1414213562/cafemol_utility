#include"FileOpenJudge.hpp"

void caf::FileOpenJudge::SuffixJudge(const std::string& filename, const std::string& suffix_name) {

	std::size_t filename_index = 0;
	std::string suffix = "." + suffix_name;

	for (const char str : filename) {

		if (str == '.') break;

		else if (filename_index < filename.size()) {
			file_prefix += str;
			++filename_index;
		}

		else {
			std::cerr << "Error: The filename '" << filename << "' is not wrong." << std::endl;
			std::exit(1);
		}
	}

	for (std::size_t idx = filename_index; idx < filename.size(); ++idx) {
		file_suffix += filename[idx];
	}

	if (file_suffix != suffix) {
		std::cerr << "Error: The format of '" << filename << "' is not '" << suffix_name << "'." << std::endl;
		std::exit(1);
	}

	else {
		std::cout << "The file '" << file_prefix << file_suffix << "' has read" << std::endl;
	}

}
