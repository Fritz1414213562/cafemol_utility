#ifndef UTIL_FUNC_HPP
#define UTIL_FUNC_HPP
#include<string>
#include<cmath>
#include<vector>

namespace cafemol {

namespace library {

std::vector<std::string> split_String(const std::string& line, const char& delim, const char& ignore_char) {

	std::vector<std::string> result;
	std::string buffer;

	for (const char& Char : line) {
		if (((Char == delim) || (Char == ignore_char)) && (buffer.empty())) continue;
		else if ((Char == delim) && (!buffer.empty())) {
			result.push_back(buffer);
			buffer.clear();
		}
		else {
			buffer.push_back(Char);
		}
	}
	if (!buffer.empty()) {
		result.push_back(buffer);
	}
	return result;
}



std::vector<int> split_String2Int(const std::string& line, const char& delim) {

	std::vector<int> result;
	std::string buffer;

	for (const char& Char : line) {
		if ((Char == delim) && (buffer.empty())) continue;
		else if ((Char == delim) && (!buffer.empty())) {
			result.push_back(std::stoi(buffer));
			buffer.clear();
		}
		else {
			buffer.push_back(Char);
		}
	}
	if (!buffer.empty()) {
		result.push_back(std::stoi(buffer));
	}
	return result;
}

std::vector<std::size_t> split_String2UnsignedLong(const std::string& line, const char& delim) {

	std::vector<std::size_t> result;
	std::string buffer;

	for (const char& Char : line) {
		if ((Char == delim) && (buffer.empty())) continue;
		else if ((Char == delim) && (!buffer.empty())) {
			result.push_back(std::stoi(buffer));
			buffer.clear();
		}
		else {
			buffer.push_back(Char);
		}
	}
	if (!buffer.empty()) {
		result.push_back(std::stoi(buffer));
	}
	return result;
}



bool is_on_same_pixel(const float& lhs_x, const float& lhs_y, const float& rhs_x, const float& rhs_y, const float& bin_width) {
	return ((((lhs_x - bin_width / 2) <= rhs_x) && (rhs_x < (lhs_x + bin_width / 2))) && 
		   (((lhs_y - bin_width / 2) <= rhs_y) && (rhs_y < (lhs_y + bin_width / 2))));
}




bool is_similar20vec(const std::vector<float>& vec, const float& cutoff) {
	bool result = true;

	for (const float& component : vec) {
		result = (result && (component <= cutoff));
	}
	return result;
}



template<typename realT>
const realT calc_EuclideanSquareDistance(const std::vector<realT>& lhs, const std::vector<realT>& rhs) {
	realT result = 0;
//	const std::size_t& iterate_num = std::min(lhs.size(), rhs.size());
	const std::size_t& iterate_num = lhs.size();

	for (std::size_t index = 0; index < iterate_num; ++index) {
		result += (lhs[index] - rhs[index]) * (lhs[index] - rhs[index]);
	}
	return result;
}


template<typename realT>
const std::vector<realT> calc_GeometricCentroid(const std::vector<std::vector<realT>>& data) {
	const std::size_t& data_number = data.size();
	const std::size_t& data_size = data[0].size();

	std::vector<realT> result(data_size, 0);
	for (std::size_t data_index = 0; data_index < data_number; ++data_index) {
		for (std::size_t index = 0; index < data_size; ++index) {
			result[index] += data[data_index][index];
		}
	}
	for (std::size_t index = 0; index < data_size; ++index) {
		result[index] /= static_cast<realT>(data_number);
	}

	return result;
}



}
}


#endif /* UTIL_FUNC_HPP */
