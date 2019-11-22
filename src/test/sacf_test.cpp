#include"../include/Statistics/StatisticsFunction.hpp"
#include<memory>
#include<array>
#include<vector>
#include<fstream>
#include<string>
#include<algorithm>
#include<cmath>

int main(int argc, char *argv[]) {

	std::string s_steps = argv[1];
	std::size_t time_steps = stoi(s_steps);

	float PI = std::acos(-1);

	std::ofstream sample_file("timedata.txt");
	std::unique_ptr<cafemol::StatisticsFunction> sacf = std::make_unique<cafemol::StatisticsFunction>();
	std::vector<float> vec(time_steps, 0);
	for (std::size_t idx = 0; idx < time_steps; ++idx) {
		float angle = idx * PI / 180;
		float size = std::log(idx + 1);
		vec[idx] = size * std::cos(angle);
		sample_file << size * std::cos(angle) << std::endl;
	}

	std::ofstream res_file("sacf_data.txt");
	std::vector<float> result = sacf->calc_SampleAutoCorrelation(vec);
	for (const float& res : result) {
		res_file << res << std::endl;
	}

	sample_file.close();
	res_file.close();

	return 0;
}
