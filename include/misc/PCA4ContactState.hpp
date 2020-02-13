#ifndef PCA_FOR_CONTACTSTATE_HPP
#define PCA_FOR_CONTACTSTATE_HPP
#include<analysis/PCA.hpp>
#include<OtherFormat/ContactStateReader.hpp>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<vector>
#include<string>
#include<array>


namespace cafemol::analysis {


class PCA4CS : public CafeInLess::analysis::PCA {

public:
	PCA4CS() : PCA() {}
	~PCA4CS() = default;


	void run(const std::string


};

}


#endif /* PCA_FOR_CONTACTSTATE_HPP */
