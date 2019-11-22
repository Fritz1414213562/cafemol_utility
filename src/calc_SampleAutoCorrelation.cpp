#include"../include/Statistics/StatisticsFunction.hpp"
#include<memory>
#include<array>
#include<vector>
#include<fstream>
#include<string>
#include<algorithm>
#include<cmath>

int main(int argc, char *argv[]) {

	std::string input_name = argv[1];
	std::string output_name = argv[2];

	std::ifstream input_file(input_name); 
	std::ofstream output_file(output_name);
	std::unique_ptr<cafemol::StatisticsFunction> sacf = std::make_unique<cafemol::StatisticsFunction>();
	std::vector<float> vec;
	std::string buf;
	while (getline(input_file, buf)) {
		int val = stoi(buf);
		vec.push_back(val);
	}

	std::vector<float> result = sacf->calc_SampleAutoCorrelation(vec);
	for (const float& res : result) {
		output_file << res << std::endl;
	}

	output_file.close();
	input_file.close();

	return 0;
}
