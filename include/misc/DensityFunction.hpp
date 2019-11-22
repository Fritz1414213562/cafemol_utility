#ifndef DENSITY_FUNCTION_HPP
#define DENSITY_FUNCTION_HPP

#include<cmath>
#include<vector>
#include<array>
#include<algorithm>

namespace cafemol {

class DensityFunction {

public:
	DensityFunction() = default;
	~DensityFunction() = default;

	void set_DataRange(const float& data_left_lim, const float& data_right_lim);
	std::array<std::vector<float>, 2> make_Density(const std::vector<float>& vec, const float bin_range);

private:
	std::array<float, 2> data_range;

};

}

#endif // DENSITY_FUNCTION_HPP
