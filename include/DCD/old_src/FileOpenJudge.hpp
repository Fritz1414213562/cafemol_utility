#include<iostream>
#include<fstream>
#include<string>

namespace caf {

class FileOpenJudge {

public:
	FileOpenJudge() = default;
	~FileOpenJudge() = default;
	std::string getFilePrefix() {return file_prefix;}
	std::string getFileSuffix() {return file_suffix;}

	void SuffixJudge(const std::string& filename, const std::string& suffix_name);

private:

	std::string file_prefix;
	std::string file_suffix;

};
}
