#include"DCDAnalyzer.hpp"

caf::DCDAnalyzer::DCDAnalyzer() : DCDParser() {
    frame_num = 0;
    atom_num = 0;
}


caf::DCDAnalyzer::DCDAnalyzer(const std::string input_file_name) : DCDParser(input_file_name) {
    frame_num = 0;
    atom_num = 0;
}


// This method can be used only when inputfile_stream was opened but not closed.
void caf::DCDAnalyzer::read_num_frame_and_atom() {

    open_file();
    // read 1st block in header and get the number of frames
    std::string first_block = read_block();
    frame_num = read_frame_num(first_block);
    // read 2nd block in header
    read_block();
    // read 3rd block in header and get the number of atoms
    std::string third_block = read_block();
    atom_num = read_atom_num(third_block);
    // option: output to standard
    std::cout << "frame_num" << frame_num << std::endl;
    std::cout << "atom_num" << atom_num << std::endl;
}


std::array<float, 3> caf::DCDAnalyzer::calcRelativeVec(const std::array<std::vector<float>, 3>& xyzs, int atom_id_begin, int atom_id_end) {

	std::array<float, 3> result;
	// take vector of begin and end
	for (std::size_t idim = 0; idim < 3; ++idim) {
		float vec_begin = xyzs[idim][atom_id_begin - 1];
		float vec_end = xyzs[idim][atom_id_end - 1];
		result[idim] = vec_end - vec_begin;
	}

	return result;
}


std::array<float, 3> caf::DCDAnalyzer::convertVec2Unit(const std::array<float, 3>& vec) {
	std::array<float, 3> result;
	std::array<float, 3> vec2 = vec * vec;
	float vec_scholar = std::sqrt(std::accumulate(vec2.begin(), vec2.end(), 0.0));
	result = vec / vec_scholar;
	return result;
}


std::array<std::array<float, 3>, 3> caf::DCDAnalyzer::adjustRotMatrix(const std::array<float, 3>& reference_vec, const std::array<float, 3>& current_vec) {

	std::array<std::array<float, 3>, 3> result;

	// adjust vector to the unit
	std::array<float, 3> ref_unit = convertVec2Unit(reference_vec);
	std::array<float, 3> cur_unit = convertVec2Unit(current_vec);

	std::array<float, 3> ref_cur_mul = ref_unit * cur_unit;
	// calculate inner_product or the angle between reference and current vector.
	float cos_bw_ref_cur = std::accumulate(ref_cur_mul.begin(), ref_cur_mul.end(), 0.0);
	float angle = std::acos(cos_bw_ref_cur);

	std::array<float, 3> normal_vec = Outer_pro(ref_unit, cur_unit);
	normal_vec = convertVec2Unit(normal_vec);
	// calculate rot
	result = RodriguesRot(normal_vec, - angle);

	return result;
}


std::vector<std::array<float, 3>> caf::DCDAnalyzer::makeVector(int atom_id_begin, int atom_id_end) {

	read_num_frame_and_atom();
	std::vector<std::array<float, 3>> result;

	for (std::size_t frame_i = 0; frame_i < frame_num; ++frame_i) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);

		std::array<float, 3> relative_vector = calcRelativeVec(xyz, atom_id_begin, atom_id_end);

		result.push_back(relative_vector);

	}

	std::cout << "relative_vector size : " << result.size() << std::endl;
	close_file();

	return result;
}


std::array<std::vector<std::array<float, 3>>, 3> caf::DCDAnalyzer::makeVector(const std::array<int, 3>& atom_id_begins, const std::array<int, 3>& atom_id_ends) {

    read_num_frame_and_atom();
    std::array<std::vector<std::array<float, 3>>, 3> result;
    
    for (std::size_t frame_i = 0; frame_i < frame_num; ++frame_i) {
        std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);

        std::array<float, 3> relative_vector = calcRelativeVec(xyz, atom_id_begins[0], atom_id_ends[0]);
		std::array<float, 3> adjust_vector = calcRelativeVec(xyz, atom_id_begins[1], atom_id_ends[1]);
		std::array<float, 3> vector_for_norm = calcRelativeVec(xyz, atom_id_begins[2], atom_id_ends[2]);
        // convert the relative vector to unit vector
		std::array<float, 3> relunit_vector = convertVec2Unit(relative_vector);
		std::array<float, 3> adjunit_vector = convertVec2Unit(adjust_vector);
		std::array<float, 3> vfnunit_vector = convertVec2Unit(vector_for_norm);

		// decomposit vfnunit_vector into adjunit_vector and norm_vector
		std::array<float, 3> adj_vfn_mul = adjunit_vector * vfnunit_vector;
		float inner_pro = std::accumulate(adj_vfn_mul.begin(), adj_vfn_mul.end(), 0.0);
		std::array<float, 3> norm_vector = vfnunit_vector - (adjunit_vector * inner_pro);
		std::array<float, 3> norunit_vector = convertVec2Unit(norm_vector);

        result[0].push_back(relunit_vector);
		result[1].push_back(adjunit_vector);
		result[2].push_back(norunit_vector);
    }

    std::cout << "relative_vector size : "  << result[0].size() << std::endl;
	std::cout << "adjust_vector size : " << result[1].size() << std::endl;
	std::cout << "normal_vector size : " << result[2].size() << std::endl;

    close_file();

    return result;
}


int caf::DCDAnalyzer::searchNativeContact(const std::array<std::vector<float>, 3>& xyzs, const std::array<int, 2>& search_chain, const std::array<int, 2>& target_chain, const float cutoff) {

	int search_begin = search_chain[0];
	int search_end = search_chain[1];
	int target_begin = target_chain[0];
	int target_end = target_chain[1];

	int result = -1;

	for (int search_index = search_begin; search_index < search_end; ++search_index) {
		for (int target_index = target_begin; target_index < target_end; ++target_index) {
			std::array<float, 3> rel_vector = calcRelativeVec(xyzs, search_index, target_index);
			std::array<float, 3> rel2 = rel_vector * rel_vector;
			float rel_scholar = std::sqrt(std::accumulate(rel2.begin(), rel2.end(), 0.0));
			if (rel_scholar < cutoff) return target_index;
			//if (rel_scholar < cutoff) return search_index;
		}
	}

	return result;
}


int caf::DCDAnalyzer::searchNativeContact(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& search_chain, const std::vector<int>& target_chain, const float cutoff) {

	int result = -1;

	for (const int& search_index : search_chain) {
		for (const int target_index : target_chain) {
			std::array<float, 3> rel_vector = calcRelativeVec(xyzs, search_index, target_index);
			std::array<float, 3> rel2 = rel_vector * rel_vector;
			float rel_scholar = std::sqrt(std::accumulate(rel2.begin(), rel2.end(), 0.0));
			if (rel_scholar < cutoff) return target_index;
			//if (rel_scholar < cutoff) return search_index;
		}
	}

	return result;
}


std::vector<int> caf::DCDAnalyzer::searchNativeContact(const std::array<int, 2>& search_chain, const std::array<int, 2>& target_chain, const float cutoff) {

	read_num_frame_and_atom();
	std::vector<int> result;

	for (std::size_t frame_i = 0; frame_i < frame_num; ++frame_i) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		int contact_index = searchNativeContact(xyz, search_chain, target_chain, cutoff);
		result.push_back(contact_index);
	}

	return result;
}


std::vector<int> caf::DCDAnalyzer::searchNativeContact(const std::array<int, 2>& search_chain, const std::array<int, 2>& target_chain, std::vector<int>& except_indices, const float cutoff) {

	read_num_frame_and_atom();
	std::vector<int> result;
	std::sort(except_indices.begin(), except_indices.end());
	int except_idx = 0;
	std::vector<int> target_ids;
	std::vector<int> search_ids;
	for (int idx = search_chain[0]; idx <= search_chain[1]; ++idx) {
		if (idx == except_indices[except_idx]) {
			++except_idx;
			continue;
		}
		search_ids.push_back(idx);
	}


	for (int idx = target_chain[0]; idx <= target_chain[1]; ++idx) {
		if (idx == except_indices[except_idx]) {
			++except_idx;
			continue;
		}
		target_ids.push_back(idx);
	}

	for (std::size_t frame_i = 0; frame_i < frame_num; ++frame_i) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		int contact_index = searchNativeContact(xyz, search_ids, target_ids, cutoff);
		result.push_back(contact_index);
	}

	return result;
}
