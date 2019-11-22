#include"../include/Statistics/GaussianFitting.hpp"
#include"../include/misc/DensityFunction.hpp"
#include<memory>
#include<array>
#include<vector>
#include<random>
#include<fstream>
#include<string>
#include<algorithm>

int main(int argc, char *argv[]) {

	std::mt19937_64 engine(10000042);

	std::normal_distribution<> dist(0.0, 10);

	std::ofstream sample_file("normal_distribution.txt");
	std::ofstream fitting_file("fitting_curve.txt");
	std::string s_max_loop = argv[1];
	std::size_t max_loop = stoi(s_max_loop);
	std::size_t steps = max_loop / 1000;

	std::vector<float> data_y(max_loop, 0.);

	for (std::size_t n = 0; n < max_loop; ++n) {
		float res = dist(engine);
		data_y[n] = res;
	}

	std::sort(data_y.begin(), data_y.end());

	std::unique_ptr<cafemol::DensityFunction> density_make = std::make_unique<cafemol::DensityFunction>();
	std::array<std::vector<float>, 2> density_data = density_make->make_Density<float>(data_y, steps);

	// std::vector<float>::iterator y_min_itr = std::min_element(data_y.begin(), data_y.end());
	// std::vector<float>::iterator y_max_itr = std::max_element(data_y.begin(), data_y.end());
	// float y_min = *y_min_itr;
	// float y_max = *y_max_itr;
	// float delta = (y_max - y_min) / steps;

	// std::vector<float> data_range;
	// std::vector<float> probability;

	// float y_point = y_min;
	// int count = 0;
	// for (std::size_t n = 0; n < max_loop; ++n) {
	// 	if ((y_point <= data_y[n]) && ((y_point + delta) > data_y[n])) {
	// 		++count;
	// 	}
	// 	else if (y_point + delta <= data_y[n]) {
	// 		float prob = static_cast<float>(count);
	// 		probability.push_back(prob);
	// 		data_range.push_back(y_point);
	// 		y_point += delta;
	// 		count = 1;
	// 	}
	// }

	for (std::size_t i = 0; i < density_data[0].size(); ++i) {
		sample_file << density_data[0][i] << " " << density_data[1][i] << std::endl;
	}

 	sample_file.close();
 
 	std::unique_ptr<cafemol::GaussianFitting> fitting = std::make_unique<cafemol::GaussianFitting>();
 	std::array<float, 3> parameters = fitting->fit_Curve(density_data[0], density_data[1]);

 	std::cout << "gaussian peak" << std::endl;
 	std::cout << parameters[0] << std::endl;
 	std::cout << "gaussian mean" << std::endl;
 	std::cout << parameters[1] << std::endl;
 	std::cout << "gaussian deviation" << std::endl;
 	std::cout << parameters[2] << std::endl;
 
 	std::vector<float> result(max_loop, 0.);
// 	for (std::size_t n = 30; n < 180; ++n) {
	for (const int n : density_data[0]) {
		// float f_n = static_cast<float>(n);
 		// float res_x = static_cast<float>(f_n / 100);
 		float res_x = static_cast<float>(n);
 		float res_y = parameters[0] 
 					* std::exp(- (res_x - parameters[1]) * (res_x - parameters[1]) 
 					/ (2 * parameters[2] * parameters[2]));
 		fitting_file << res_x << " " << res_y << std::endl;
 	}
 
 	fitting_file.close();

	return 0;
}
