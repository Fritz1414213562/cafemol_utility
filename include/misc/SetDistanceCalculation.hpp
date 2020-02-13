#ifndef SET_DISTANCE_CALCULATION_HPP
#define SET_DISTANCE_CALCULATION_HPP
#include<vector>
#include<opencv2/core.hpp>
#include<opencv2/imgproc.hpp>



namespace cafemol::library {


const float calc_EarthMoversDistance(const std::vector<cv::Vec3f>& lhs, const std::vector<cv::Vec3f>& rhs) {
	const cv::Mat& lhs_mat = cv::Mat(lhs).clone().reshape(1);
	const cv::Mat& rhs_mat = cv::Mat(rhs).clone().reshape(1);
	float result = cv::EMD(lhs_mat, rhs_mat, cv::DIST_L2);

	return result;
}



}





#endif /* SET_DISTANCE_CALCULATION_HPP */
