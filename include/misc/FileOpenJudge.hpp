#include<iostream>
#include<fstream>
#include<array>
#include<string>

namespace cafemol {

class FileOpenJudge {

public:
	FileOpenJudge() = default;
	~FileOpenJudge() = default;
	std::string getFilePrefix() {return file_prefix;}
	std::string getFileSuffix() {return file_suffix;}

	void SuffixJudge(const std::string& filename, const std::string& suffix_name);

//	template<std::size_t file_num>
//	std::array<std::string, file_num> SuffixesJudge(const std::array<std::string, file_num>& filenames, const std::array<std::string, file_num>& suffix_names);

	template<std::size_t file_num>
	std::array<std::string, file_num> SuffixesJudge(const std::array<std::string, file_num>& filenames, const std::array<std::string, file_num>& suffix_names) {
	
		std::array<std::string, file_num> result;
	
		std::array<std::string, file_num> answer_suffixes;
		for (std::size_t idx = 0; idx < file_num; ++idx) {
			answer_suffixes[idx] = "." + suffix_names[idx];
		}
	
		std::array<std::string, file_num> query_suffixes;
		for (std::size_t idx = 0; idx < file_num; ++idx) {
			const std::string& filename = filenames[idx];
			std::string reverse_suffix;
			std::string suffix;
			std::size_t filename_length = filename.size();
			int reverse_idx = filename_length;
	
			// search the suffix of filename
			for (int filename_idx = filename_length - 1; filename_idx >= 0; --filename_idx) {
				const char suffix_str = filename[filename_idx];
	
				if (suffix_str == '.') {
					reverse_suffix.push_back(suffix_str);
					--reverse_idx;
					break;
				}
				else if (idx < 0) {
					std::cerr << "Error: Define the file format." << std::endl;
					std::cerr << filename << " <-" << std::endl;
					std::exit(1);
				}
				else {
					reverse_suffix.push_back(suffix_str);
					--reverse_idx;
				}
			}
	
			// reverse reverse_suffix
			for (int rev_idx = reverse_suffix.size() - 1; rev_idx >= 0; --rev_idx) {
				suffix.push_back(reverse_suffix[rev_idx]);
			}
	
			// search the suffix consistent with filesuffix
			std::size_t result_idx = 0;
			for (const std::string& answer_suffix : answer_suffixes) {
				if (answer_suffix != suffix) {
					++result_idx;
					continue;
				}
				else break;
			}
	
			if (result_idx >= file_num) continue; //{
	//			std::cerr << "There is a strange format file." << std::endl;
	//			std::cerr << "The file is not formatted as these formats." << std::endl;
	//			for (const std::string& suffix_name : suffix_names) {
	//				std::cerr << suffix_name << std::endl;
	//			}
	//			std::exit(1);
	//		}
			else if (!result[result_idx].empty()) {
				std::cerr << "Something wrong! The index is overlapped." << std::endl;
				std::cerr << "Probably, two file exists specified as the same format." << std::endl;
				std::exit(1);
			}
			else {
				result[result_idx] = filename;
			}
		}
	
		int check_idx = 0;
		for (const std::string& result_filename : result) {
			if (result_filename.empty()) {
				std::cerr << "Error: The strange file exists." << std::endl;
				std::cerr << "Please use this file format." << std::endl;
				std::cerr << " -> " << suffix_names[check_idx] << std::endl;
				std::exit(1);
			}
			++check_idx;
		}
	
		return result;
	}

	// operator overload for objective function
	void operator()(const std::string& filename, const std::string& suffix_name) {
		SuffixJudge(filename, suffix_name);
	}

	template<std::size_t file_num>
	std::array<std::string, file_num> operator()(const std::array<std::string, file_num>& filenames, const std::array<std::string, file_num>& suffix_names) {
		return SuffixesJudge(filenames, suffix_names);
	}

private:

	std::string file_prefix;
	std::string file_suffix;

	void clear_Prefix() {if (!file_prefix.empty()) file_prefix.clear();}
	void clear_Suffix() {if (!file_suffix.empty()) file_suffix.clear();}

};
}
