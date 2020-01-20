#ifndef DETECT_CONTACTS_HPP
#define DETECT_CONTACTS_HPP
#include<CafemolName.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<iostream>
#include<string>
#include<array>
#include<vector>


namespace cafemol {

class ContactSearcher {

public:
	ContactSearcher() = default;
	ContactSearcher(const float& cutoff) : cutoff_length(cutoff) {}
	~ContactSearcher() = default;

	void set_CutoffLength(const float& cutoff);

	std::vector<int> get_TimeSeriesOfClosestID2(const std::size_t& atom_id, const Trajectory& traj);
	std::vector<int> get_TimeSeriesOfClosestID2(const std::size_t& atom_id, const std::vector<std::size_t>& iterate_range, const Trajectory& traj);

	std::vector<std::vector<float>> get_DistanceMatrix(const std::array<std::vector<float>, 3>& snapshot);

	std::array<std::array<int, 2>, 2> get_ContactEnds(const std::vector<std::size_t>& source_range, const std::vector<std::size_t>& target_range, const std::array<std::vector<float>, 3>& snapshot);


private:
	// default
	float cutoff_length = 10.0;

	// member object
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();

	int get_ClosestID2(const std::size_t& atom_id, const std::array<std::vector<float>, 3>& snapshot);
	int get_ClosestID2(const std::size_t& atom_id, const std::vector<std::size_t>& iterate_range, const std::array<std::vector<float>, 3>& snapshot);

	std::array<float, 3> get_Vector(const std::size_t& atom_index, const std::array<std::vector<float>, 3>& snapshot);
};
}


#endif /* DETECT_CONTACTS_HPP */
