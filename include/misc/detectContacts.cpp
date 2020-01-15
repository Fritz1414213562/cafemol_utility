#include<detectContacts.hpp>


void cafemol::ContactSearcher::set_CutoffLength(const float& cutoff) {
	cutoff_length = cutoff;
}


std::vector<int> cafemol::ContactSearcher::get_TimeSeriesOfClosestID2(const std::size_t& atom_id, const Trajectory& traj) {

	std::vector<int> result;

	for (const std::array<std::vector<float>, 3>& snapshot : traj) {
		const int& closest_id = get_ClosestID2(atom_id, snapshot) + 1;
		result.push_back(closest_id);
	}
	return result;
}


std::vector<int> cafemol::ContactSearcher::get_TimeSeriesOfClosestID2(const std::size_t& atom_id, const std::vector<std::size_t>& iterate_range, const Trajectory& traj) {

	std::vector<int> result;

	for (const std::array<std::vector<float>, 3>& snapshot : traj) {
		const int& closest_id = get_ClosestID2(atom_id, iterate_range, snapshot) + 1;
		result.push_back(closest_id);
	}
	return result;
}


std::vector<std::vector<float>> cafemol::ContactSearcher::get_DistanceMatrix(const std::array<std::vector<float>, 3>& snapshot) {

	if ((snapshot[0].size() != snapshot[1].size()) || (snapshot[0].size() != snapshot[2].size())) {
		eout("The degree of freedom at each direction is not consistent. This trajectory might be broken.");
	}
	const std::size_t& atom_size = snapshot[0].size();

	std::vector<std::vector<float>> result(atom_size);

	for (std::size_t i_atom = 0; i_atom < atom_size - 1; ++i_atom) {

		result[i_atom].resize(atom_size, -1.0);
		const std::array<float, 3>& vector_index_i = get_Vector(i_atom, snapshot);

		for (std::size_t j_atom = i_atom; j_atom < atom_size; ++j_atom) {
			const std::array<float, 3>& vector_index_j = get_Vector(j_atom, snapshot);
			// calculate a distance between i and j
			float dist2 = 0.0;
			for (std::size_t idim = 0; idim < 3; ++idim) {
				dist2 += (vector_index_i[idim] - vector_index_j[idim]) * 
						 (vector_index_i[idim] - vector_index_j[idim]);
			}
			if (dist2 > cutoff_length * cutoff_length) continue;
			else result[i_atom][j_atom] = dist2;
		}
	}
	return result;

}


int cafemol::ContactSearcher::get_ClosestID2(const std::size_t& atom_id, const std::array<std::vector<float>, 3>& snapshot) {

	if ((snapshot[0].size() != snapshot[1].size()) || (snapshot[0].size() != snapshot[2].size())) {
		eout("The degree of freedom at each direction is not consistent. This trajectory might be broken.");
	}
	const std::size_t& atom_size = snapshot[0].size();

	if ((1 > atom_id) || (atom_id > atom_size)) eout("The AtomID is out of range");

	const std::array<float, 3>& source_atom_vector = get_Vector(atom_id - 1, snapshot);

	float min_dist2 = cutoff_length * cutoff_length;
	int min_idx = -1;

	for (std::size_t i_atom = 0; i_atom < atom_size; ++i_atom) {
		// cycle when comparing between the same ID.
		if (i_atom == atom_id - 1) continue;

		const std::array<float, 3>& target_atom_vector = get_Vector(i_atom, snapshot);
		float dist2 = 0.0;
		for (std::size_t idim = 0; idim < 3; ++idim) {
			dist2 += (target_atom_vector[idim] - source_atom_vector[idim]) * (target_atom_vector[idim] - source_atom_vector[idim]);
		}
		if (dist2 > cutoff_length * cutoff_length) continue;
		else if (dist2 <= min_dist2) {
			min_dist2 = dist2;
			min_idx = i_atom;
		}
	}
	return min_idx;
}


int cafemol::ContactSearcher::get_ClosestID2(const std::size_t& atom_id, const std::vector<std::size_t>& iterate_range, const std::array<std::vector<float>, 3>& snapshot) {

	if ((snapshot[0].size() != snapshot[1].size()) || (snapshot[0].size() != snapshot[2].size())) {
		eout("The degree of freedom at each direction is not consistent. This trajectory might be broken.");
	}
	const std::size_t& atom_size = snapshot[0].size();

	if ((1 > atom_id) || (atom_id > atom_size)) eout("The AtomID is out of range");

	const std::array<float, 3>& source_atom_vector = get_Vector(atom_id - 1, snapshot);

	float min_dist2 = cutoff_length * cutoff_length;
	int min_idx = -1;

	for (const std::size_t& i_atom : iterate_range) {
		// cycle when comparing between the same ID.
		if (i_atom == atom_id) continue;
		else if ((1 > i_atom) || (i_atom > atom_size)) eout("The iterate range is out of atom id range");

		const std::array<float, 3>& target_atom_vector = get_Vector(i_atom - 1, snapshot);
		float dist2 = 0.0;
		for (std::size_t idim = 0; idim < 3; ++idim) {
			dist2 += (target_atom_vector[idim] - source_atom_vector[idim]) * (target_atom_vector[idim] - source_atom_vector[idim]);
		}
		if (dist2 > cutoff_length * cutoff_length) continue;
		else if (dist2 <= min_dist2) {
			min_dist2 = dist2;
			min_idx = i_atom - 1;
		}
	}
	return min_idx;
}


std::array<float, 3> cafemol::ContactSearcher::get_Vector(const std::size_t& atom_index, const std::array<std::vector<float>, 3>& snapshot) {
	const std::array<float, 3>& res = {snapshot[0][atom_index],
									   snapshot[1][atom_index],
									   snapshot[2][atom_index]};
	return res;
}
