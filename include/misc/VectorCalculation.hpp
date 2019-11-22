#ifndef VECTOR_CALCULATION_HPP
#define VECTOR_CALCULATION_HPP

#include<array>
#include<cmath>


namespace cafemol {

namespace library {


template<typename realT, std::size_t sizeN>
realT calc_Distance(const std::array<realT, sizeN>& lhs, const std::array<realT, sizeN>& rhs) {
	realT result_sqr = 0;
	for (std::size_t idx = 0; idx < sizeN; ++idx) {
		result_sqr += (lhs[idx] - rhs[idx]) * (lhs[idx] - rhs[idx]);
	}
	realT result = std::sqrt(result_sqr);
	return result;
}

}
}


#endif /* VECTOR_CALCULATION_HPP */
