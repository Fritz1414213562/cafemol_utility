#ifndef CAFEMOL_NAME_HPP
#define CAFEMOL_NAME_HPP
#include<Eigen/Core>

namespace cafemol {

// common alias of cafemol
using Trajectory = std::vector<std::array<std::vector<float>, 3>>;
using EigenTrajectory = std::vector<Eigen::VectorXd>;

namespace library {
// alias for library namespace
}

}

#endif /* CAFEMOL_NAME_HPP */
