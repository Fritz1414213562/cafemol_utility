#ifndef PARSECOMMANDLINEARGUMENTS_HPP
#define PARSECOMMANDLINEARGUMENTS_HPP
#include<vector>
#include<string>


namespace cafemol {


namespace library {


std::vector<std::string> split_ByHyphen(const std::string& input_string) {

	std::vector<std::string> output_strings;

	std::string buffer;
	for (const char& Char : input_string) {
		if (Char == '-') {
			output_strings.push_back(buffer);
			buffer.clear();
			continue;
		}
		buffer.push_back(Char);
	}
	if (buffer.empty()) {
		std::cerr << "The input string is wrong. Probably, '-' comes to the ends of input string." << std::endl;
		std::exit(1);
	}

	output_strings.push_back(buffer);

	return output_strings;
}


}
}


#endif /* PARSECOMMANDLINEARGUMENTS_HPP */
