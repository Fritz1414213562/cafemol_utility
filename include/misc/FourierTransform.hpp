#ifndef FOURIER_TRANSFORM_HPP
#define FOURIER_TRANSFORM_HPP

#include<array>
#include<vector>
#include<cmath>
#include<iostream>
#include<string>


namespace cafemol {

class FourierTransform {

public:
	FourierTransform() = default;
	~FourierTransform() = default;

protected:
	
	template<typename realT>
	std::vector<realT> calc_DisperseFourierTransform(const std::vector<realT>& vec) {
		for (std::size_t idx = 0; ) {
			for (std::size_t freq = 0;
		}
	}

private:
	const float frequency_max = 2 * std::acos(-1);

};

}

#endif /* FOURIER_TRANSFORM_HPP */
