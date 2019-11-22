#include"GaussianFitting.hpp"


void cafemol::GaussianFitting::initialize_FittingMatrix(const std::vector<float>& data_x, const std::vector<float>& data_y) {
	
	// calculate the components of fitting matrix.
	float y_sqr_sum;
	float x_y_sqr_sum;
	float x_sqr_y_sqr_sum;
	float x_cubic_y_sqr_sum;
	float x_forth_y_sqr_sum;

	{
		std::vector<float> y_sqr = data_y * data_y;
		y_sqr_sum = std::accumulate(y_sqr.begin(), y_sqr.end(), 0.0);

		std::vector<float> x_y_sqr = data_x * y_sqr;
		x_y_sqr_sum = std::accumulate(x_y_sqr.begin(), x_y_sqr.end(), 0.0);

		std::vector<float> x_sqr_y_sqr = data_x * x_y_sqr;
		x_sqr_y_sqr_sum = std::accumulate(x_sqr_y_sqr.begin(), x_sqr_y_sqr.end(), 0.0);

		std::vector<float> x_cubic_y_sqr = data_x * x_sqr_y_sqr;
		x_cubic_y_sqr_sum = std::accumulate(x_cubic_y_sqr.begin(), x_cubic_y_sqr.end(), 0.0);

		std::vector<float> x_forth_y_sqr = data_x * x_cubic_y_sqr;
		x_forth_y_sqr_sum = std::accumulate(x_forth_y_sqr.begin(), x_forth_y_sqr.end(), 0.0);
	}

	fitting_matrix << y_sqr_sum,		x_y_sqr_sum,		x_sqr_y_sqr_sum, 
					  x_y_sqr_sum,		x_sqr_y_sqr_sum,	x_cubic_y_sqr_sum,
					  x_sqr_y_sqr_sum,	x_cubic_y_sqr_sum,	x_forth_y_sqr_sum; 

}


void cafemol::GaussianFitting::initialize_OutputsVector(const std::vector<float>& data_x, const std::vector<float>& data_y) {

	// calculate the components of outputs.
	float outputs_x;
	float outputs_y;
	float outputs_z;

	{
		std::vector<float> y_sqr = data_y * data_y;
		std::vector<float> log_y(data_y.size(), 0);
		for (std::size_t idx = 0; idx < data_y.size(); ++idx) {
			log_y[idx] = std::log(data_y[idx]);
		}
		std::vector<float> y_sqr_log_y = y_sqr * log_y;
		std::vector<float> x_y_sqr_log_y = data_x * y_sqr_log_y;
		std::vector<float> x_sqr_y_sqr_log_y = data_x * x_y_sqr_log_y;

		outputs_x = std::accumulate(y_sqr_log_y.begin(), y_sqr_log_y.end(), 0.0);
		outputs_y = std::accumulate(x_y_sqr_log_y.begin(), x_y_sqr_log_y.end(), 0.0);
		outputs_z = std::accumulate(x_sqr_y_sqr_log_y.begin(), x_sqr_y_sqr_log_y.end(), 0.0);
	}

	outputs_vector << outputs_x, outputs_y, outputs_z;

}


Eigen::Vector3f cafemol::GaussianFitting::solve_LinearSystem() {

	Eigen::Vector3f result = fitting_matrix.fullPivLu().solve(outputs_vector);
	return result;
}


std::array<float, 3> cafemol::GaussianFitting::fit_Curve(const std::vector<float>& data_x, const std::vector<float>& data_y) {

	// initialize fitting matrix and outputs vector
	initialize_FittingMatrix(data_x, data_y);
	initialize_OutputsVector(data_x, data_y);

	Eigen::Vector3f inputs_vector = solve_LinearSystem();
	float peak_Gaussian = std::exp(inputs_vector(0)
								 - inputs_vector(1) * inputs_vector(1) / (4 * inputs_vector(2)));

	float mean_Gaussian = - inputs_vector(1) / (2 * inputs_vector(2));
	float deviation_Gaussian = std::sqrt(- 1 / (2 * inputs_vector(2)));

	std::array<float, 3> result = {peak_Gaussian, mean_Gaussian, deviation_Gaussian};
	return result;
}
