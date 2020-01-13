#include"DCDAnalyzer.hpp"

cafemol::DCDAnalyzer::DCDAnalyzer() : DCDParser() {
    frame_num = 0;
    atom_num = 0;
}


cafemol::DCDAnalyzer::DCDAnalyzer(const std::string input_file_name) : DCDParser(input_file_name) {
    frame_num = 0;
    atom_num = 0;
}


int cafemol::DCDAnalyzer::get_AtomNumber() {
	return atom_num;
}


std::array<float, 3> cafemol::DCDAnalyzer::getVector(const std::array<std::vector<float>, 3>& xyzs, const int vec_id) {

	std::array<float, 3> result;
	// take vector at the id
	for (std::size_t idim = 0; idim < 3; ++idim) {
		float component = xyzs[idim][vec_id - 1];
		result[idim] = component;
	}

	return result;
}


std::array<float, 3> cafemol::DCDAnalyzer::convertVec2Unit(const std::array<float, 3>& vec) {
	std::array<float, 3> result;
	std::array<float, 3> vec2 = vec * vec;
	float vec_scholar = std::sqrt(std::accumulate(vec2.begin(), vec2.end(), 0.0));
	result = vec / vec_scholar;
	return result;
}


std::array<float, 3> cafemol::DCDAnalyzer::calcCenter(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& ids) {
	std::array<float, 3> result = {0., 0., 0.};

	int ids_num = ids.size();
	for (const std::size_t& id : ids) {
		std::array<float, 3> vector_at_id = getVector(xyzs, id);
		result += vector_at_id; 
	}
	float ids_num_real = static_cast<float>(ids_num);
	result /= ids_num_real;

	return result;
}


int cafemol::DCDAnalyzer::getClosestResidue(const std::array<std::vector<float>, 3>& xyzs, const std::array<float, 3>& fixed_points_vector, const std::vector<int>& search_chain_ids) {

	// initial value
	int result = -1;
	std::array<float, 3> ini_vector = getVector(xyzs, search_chain_ids[0]);
	float min_distance = calcDistance(fixed_points_vector, ini_vector);
	// search the closest residue
	for (const int& search_id : search_chain_ids) {
		std::array<float, 3> vector_at_id = getVector(xyzs, search_id);
		float distance = calcDistance(fixed_points_vector, vector_at_id);
		if (distance < min_distance) {
			min_distance = distance;
			result = search_id;
		}
	}

	return result;
}


int cafemol::DCDAnalyzer::getClosestResidue(const std::array<std::vector<float>, 3>& xyzs, const std::array<float, 3>& fixed_points_vector, const std::vector<int>& search_chain_ids, const float& cutoff) {

	// initial value
	int result = -1;
	std::array<float, 3> ini_vector = getVector(xyzs, search_chain_ids[0]);
	float min_distance = calcDistance(fixed_points_vector, ini_vector);
	// search the closest residue
	for (const int& search_id : search_chain_ids) {
		std::array<float, 3> vector_at_id = getVector(xyzs, search_id);
		float distance = calcDistance(fixed_points_vector, vector_at_id);
		if (distance > cutoff) continue;
		else if (distance < min_distance) {
			min_distance = distance;
			result = search_id;
		}
	}

	return result;
}


std::vector<int> cafemol::DCDAnalyzer::getClosestResidue(const std::vector<int>& fixed_points_ids, const std::vector<int>& search_chain_ids) {

	// calculate the center of fixed points
	read_num_frame_and_atom();
	std::vector<int> result;
	result.resize(frame_num);

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		std::array<float, 3> fixed_points_center = calcCenter(xyz, fixed_points_ids);

		// for (const float& comp : fixed_points_center) {
		// 	std::cout << comp << " ";
		// }
		// std::cout << std::endl;

		int closest_id = getClosestResidue(xyz, fixed_points_center, search_chain_ids);
		result[iframe] = closest_id;
	}

	load_file();

	return result;
}


std::vector<int> cafemol::DCDAnalyzer::getClosestResidue(const std::vector<int>& fixed_points_ids, const std::vector<int>& search_chain_ids, const float& cutoff) {

	// calculate the center of fixed points
	read_num_frame_and_atom();
	std::vector<int> result;
	result.resize(frame_num);

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		std::array<float, 3> fixed_points_center = calcCenter(xyz, fixed_points_ids);

		int closest_id = getClosestResidue(xyz, fixed_points_center, search_chain_ids, cutoff);
		result[iframe] = closest_id;
	}

	load_file();

	return result;
}


std::vector<std::array<int, 2>> cafemol::DCDAnalyzer::getClosestResidue(const std::array<int, 2>& fixed_points_ids, const std::vector<int>& search_chain_ids, const float& cutoff) {

	// calculate the center of fixed points
	read_num_frame_and_atom();
	std::vector<std::array<int, 2>> result;
	result.resize(frame_num);

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		std::array<float, 3> fixed_point1_vector = {xyz[0][fixed_points_ids[0] - 1],
													xyz[1][fixed_points_ids[0] - 1],
													xyz[2][fixed_points_ids[0] - 1]};
		std::array<float, 3> fixed_point2_vector = {xyz[0][fixed_points_ids[1] - 1],
													xyz[1][fixed_points_ids[1] - 1],
													xyz[2][fixed_points_ids[1] - 1]};
		int closest_id1 = getClosestResidue(xyz, fixed_point1_vector, search_chain_ids, cutoff); 
		int closest_id2 = getClosestResidue(xyz, fixed_point2_vector, search_chain_ids, cutoff); 
		result[iframe] = {closest_id1, closest_id2};
	}

	load_file();

	return result;
}


bool cafemol::DCDAnalyzer::is_NativeContact(const std::array<float, 3>& vec1, const std::array<float, 3>& vec2, const float& cutoff) {
	
	float distance = calcDistance(vec1, vec2);
	bool result = (distance <= cutoff);
	return result;
}


std::array<int, 2> cafemol::DCDAnalyzer::get_1stNativeContactedResidue(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff) {

	std::array<int, 2> result = {-1, -1};
	bool is_searched = false;
	for (const int& search_chain_id : search_chain_ids) {
		std::array<float, 3> search_vec = {xyzs[0][search_chain_id - 1], 
										   xyzs[1][search_chain_id - 1], 
										   xyzs[2][search_chain_id - 1]};
	
		for (const int& contact_target : contact_targets) {
			std::array<float, 3> target_vec = {xyzs[0][contact_target - 1],
											   xyzs[1][contact_target - 1],
											   xyzs[2][contact_target - 1]};
			if (is_NativeContact(search_vec, target_vec, cutoff)) {
				result[0] = search_chain_id;
				result[1] = contact_target;
				is_searched = true;
				break;
			}
		}
		if (is_searched) break;
	}
	return result;
}


std::array<int, 2> cafemol::DCDAnalyzer::get_lastNativeContactedResidue(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff) {

	std::array<int, 2> result = {-1, -1};
	bool is_searched = false;
	for (int idx = search_chain_ids.size() - 1; idx >= 0; --idx) {

		int search_chain_id = search_chain_ids[idx];

		std::array<float, 3> search_vec = {xyzs[0][search_chain_id - 1],
										   xyzs[1][search_chain_id - 1],
										   xyzs[2][search_chain_id - 1]};

		for (const int& contact_target : contact_targets) {
			std::array<float, 3> target_vec = {xyzs[0][contact_target - 1],
											   xyzs[1][contact_target - 1],
											   xyzs[2][contact_target - 1]};
			if (is_NativeContact(search_vec, target_vec, cutoff)) {
				result[0] = search_chain_id;
				result[1] = contact_target;
				is_searched = true;
				break;
			}
		}
		if (is_searched) break;
	}
	return result;
}


std::vector<std::array<int, 2>> cafemol::DCDAnalyzer::get_1stNativeContactedResidue(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff) {

	read_num_frame_and_atom();
	std::vector<std::array<int, 2>> result;
	result.resize(frame_num, {0, 0});

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		std::array<int, 2> first_native_contact_id = get_1stNativeContactedResidue(xyz, search_chain_ids, contact_targets, cutoff);
		result[iframe] = first_native_contact_id;
	}

	load_file();

	return result;
}


std::array<std::vector<int>, 2> cafemol::DCDAnalyzer::get_NativeContactedRange(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff) {
	
	read_num_frame_and_atom();
	std::array<std::vector<int>, 2> result;
	result[0].resize(frame_num, 0);
	result[1].resize(frame_num, 0);

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		std::array<int, 2> first_native_contact_id = get_1stNativeContactedResidue(xyz, search_chain_ids, contact_targets, cutoff);
		if (first_native_contact_id[0] <= 0) {
			result[0][iframe] = first_native_contact_id[0];
			result[1][iframe] = -1;
		}
		else {
			std::array<int, 2> last_native_contact_id = get_lastNativeContactedResidue(xyz, search_chain_ids, contact_targets, cutoff);
			result[0][iframe] = first_native_contact_id[0];
			result[1][iframe] = last_native_contact_id[0];
		}
	}

	load_file();

	return result;
}


int cafemol::DCDAnalyzer::count_NativeContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& ini_contact_pairs, const float& cutoff) {

	int result = 0;
	for (const std::array<int, 2>& ini_contact_pair : ini_contact_pairs) {
		int mass_point1_id = ini_contact_pair[0];
		int mass_point2_id = ini_contact_pair[1];
		std::array<float, 3> mass_point1 = {xyzs[0][mass_point1_id - 1],
											xyzs[1][mass_point1_id - 1],
											xyzs[2][mass_point1_id - 1]};
		std::array<float, 3> mass_point2 = {xyzs[0][mass_point2_id - 1],
											xyzs[1][mass_point2_id - 1],
											xyzs[2][mass_point2_id - 1]};
		if (is_NativeContact(mass_point1, mass_point2, cutoff)) ++result;
	}
    return result;
}


int cafemol::DCDAnalyzer::count_NativeContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& ini_contact_pairs, const std::vector<float>& cutoff) {

    int result = 0;
    int idx_cutoff = 0;
    for (const std::array<int, 2>& ini_contact_pair : ini_contact_pairs) {
        int mass_point1_id = ini_contact_pair[0];
        int mass_point2_id = ini_contact_pair[1];
        std::array<float, 3> mass_point1 = {xyzs[0][mass_point1_id - 1],
                                            xyzs[1][mass_point1_id - 1],
                                            xyzs[2][mass_point1_id - 1]};
        std::array<float, 3> mass_point2 = {xyzs[0][mass_point2_id - 1],
                                            xyzs[1][mass_point2_id - 1],
                                            xyzs[2][mass_point2_id - 1]};
        if (is_NativeContact(mass_point1, mass_point2, cutoff[idx_cutoff])) {
            ++result;
        }
        ++idx_cutoff;
    }
    return result;

}


std::vector<int> cafemol::DCDAnalyzer::count_NativeContacts(const std::vector<std::array<int, 2>>& ini_contact_pairs, const float& cutoff) {

	read_num_frame_and_atom();
	std::vector<int> result(frame_num, 0);

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		int contact_num = count_NativeContacts(xyz, ini_contact_pairs, cutoff);
		result[iframe] = contact_num;
	}

	load_file();

	return result;
}


std::vector<int> cafemol::DCDAnalyzer::count_NativeContacts(const std::vector<std::array<int, 2>>& ini_contact_pairs, const std::vector<float>& cutoff) {

    if (ini_contact_pairs.size() != cutoff.size()) {
        std::cerr << "Error: The number of initial contact pairs is not consistent with the number of cutoffs." << std::endl;
        std::exit(1);
    }

    read_num_frame_and_atom();
    std::vector<int> result(frame_num, 0);

    for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
        std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
        int contact_num = count_NativeContacts(xyz, ini_contact_pairs, cutoff);
        result[iframe] = contact_num;
    }

	load_file();

    return result;
}



std::vector<int> cafemol::DCDAnalyzer::get_NativeContacts(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff, const int& frame) {

	read_num_frame_and_atom();

	if (frame_num < frame) {
		std::cerr << "Error: This trajectory doesn't have the snapshot of the frame '" << frame << "'" << std::endl;
		std::exit(1);
	}
	else if (frame < 1) {
		std::cerr << "Error: Frame number is smaller than 1" << std::endl;
		std::exit(1);
	}

	for (std::size_t iframe = 0; iframe < frame - 1; ++iframe) {
		read_xyz(atom_num);
	}

	std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
	std::vector<int> result;

	for (const int& search_chain_id : search_chain_ids) {
		float min_dist = 1000.0;
		int min_id = -1000;
		for (const int& contact_target : contact_targets) {
			std::array<float, 3> source_vect = {xyz[0][search_chain_id - 1],
												xyz[1][search_chain_id - 1],
												xyz[2][search_chain_id - 1]};
			std::array<float, 3> target_vect = {xyz[0][contact_target - 1],
												xyz[1][contact_target - 1],
												xyz[2][contact_target - 1]};
			float dist_source_target = calcDistance(source_vect, target_vect);
			if (dist_source_target > cutoff) continue;

			else if (min_dist > dist_source_target && !cafemol::library::is_contains<int>(result, contact_target)) {
				min_id = contact_target;
				min_dist = dist_source_target;
			}
		}
		if (min_id < 0) continue;
		result.push_back(min_id);
	}

	load_file();

	return result;
}


std::vector<int> cafemol::DCDAnalyzer::get_NativeContactsFromIni(const std::vector<int>& search_chain_ids, const std::vector<int>& contact_targets, const float& cutoff) {
	return get_NativeContacts(search_chain_ids, contact_targets, cutoff, 1);
}


std::vector<int> cafemol::DCDAnalyzer::get_AtomIDsNearPoint(const std::array<float, 3>& fixed_point_coordinate, const std::vector<int>& search_chain_ids, const float& cutoff) {

	read_num_frame_and_atom();
	std::vector<int> result(frame_num);
	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		result[iframe] = getClosestResidue(xyz, fixed_point_coordinate, search_chain_ids, cutoff);
	}
	load_file();

	return result;
}


std::vector<int> cafemol::DCDAnalyzer::count_PDNSContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float, 3>>& parameters) {
	std::vector<int> result;
	
	const std::size_t pdns_size = pdns_pro_ids.size();
	for (std::size_t idx = 0; idx < pdns_size; ++idx) {
		int pdns_ca_id = pdns_pro_ids[idx];
		int pdns_ca_N;
		int pdns_ca_C;

		std::array<float, 3> parameter = parameters[idx];

		if (pdns_ca_id < 1) error_output("Invalid Atom ID", pdns_ca_id, "negative value error");
		else if (pdns_ca_id == 1) {
			pdns_ca_N = 1;
			pdns_ca_C = pdns_ca_id + 1;
		}
		else if (pdns_ca_id == xyzs[0].size()) {
			pdns_ca_N = pdns_ca_id - 1;
			pdns_ca_C = xyzs[0].size();
		}
		else {
			pdns_ca_N = pdns_ca_id - 1;
			pdns_ca_C = pdns_ca_id + 1;
		}

		std::array<float, 3> pdns_ca_C2N_vec = getVector(xyzs, pdns_ca_N) - getVector(xyzs, pdns_ca_C);
		float pdns_ca_C2N_scholar = calc_Norm(pdns_ca_C2N_vec);

		std::array<float, 3> pdns_ca_vec = getVector(xyzs, pdns_ca_id);

		for (const std::array<int, 2>& phos_sug_id : pdns_phos_sug_ids) {
			std::array<float, 3> pdns_phos_vec = getVector(xyzs, phos_sug_id[0]);
			if (calcDistance(pdns_phos_vec, pdns_ca_vec) > (5.0 + parameter[0])) continue;

			std::array<float, 3> pdns_phos2ca_vec = pdns_ca_vec - pdns_phos_vec;
			float pdns_phos2ca_scholar = calc_Norm(pdns_phos2ca_vec);
			std::array<float, 3> pdns_phos2sug_vec = getVector(xyzs, phos_sug_id[1]) - pdns_phos_vec;
			float pdns_phos2sug_scholar = calc_Norm(pdns_phos2sug_vec);

			float cos_ca2phos2sug = calcDot(pdns_phos2ca_vec, pdns_phos2sug_vec) / (pdns_phos2ca_scholar * pdns_phos2sug_scholar);
			float cos_N2C2phos = calcDot(pdns_ca_C2N_vec, pdns_phos2ca_vec) / (pdns_ca_C2N_scholar * pdns_phos2ca_scholar);
			float rad_ca2phos2sug = acos(cos_ca2phos2sug);
			float rad_N2C2phos = acos(cos_N2C2phos);

			float d_rad_ca2phos2sug = rad_ca2phos2sug - parameter[1];
			float d_rad_N2C2phos = rad_N2C2phos - parameter[2];
			if (((-2 * cutoff_angle) <= d_rad_ca2phos2sug) && 
				(d_rad_ca2phos2sug <= 2 * cutoff_angle) &&
				((-2 * cutoff_angle) <= d_rad_N2C2phos) &&
				(d_rad_N2C2phos <= 2 * cutoff_angle)) {
				result.push_back(pdns_ca_id);
				break;
			}
		}
	}
	return result;
}


std::vector<std::vector<int>> cafemol::DCDAnalyzer::count_PDNSContacts(const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float, 3>>& parameters) {

	read_num_frame_and_atom();
	std::vector<std::vector<int>> result(frame_num);

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		result[iframe] = count_PDNSContacts(xyz, pdns_phos_sug_ids, pdns_pro_ids, parameters);
	}
	load_file();

	return result;
}


std::vector<std::array<int, 2>> cafemol::DCDAnalyzer::get_PDNSContactResidues(const std::array<std::vector<float>, 3>& xyzs, const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float ,3>>& parameters) {

	std::vector<std::array<int, 2>> result;
	std::size_t pdns_size = pdns_pro_ids.size();
	for (std::size_t idx = 0; idx < pdns_size; ++idx) {
		int pdns_ca_id = pdns_pro_ids[idx];
		int pdns_ca_N;
		int pdns_ca_C;

		std::array<float, 3> parameter = parameters[idx];

		if (pdns_ca_id < 1) error_output("Invalid Atom ID", pdns_ca_id, "negative value error");
		else if (pdns_ca_id == 1) {
			pdns_ca_N = 1;
			pdns_ca_C = pdns_ca_id + 1;
		}
		else if (pdns_ca_id == xyzs[0].size()) {
			pdns_ca_N = pdns_ca_id - 1;
			pdns_ca_C = xyzs[0].size();
		}
		else {
			pdns_ca_N = pdns_ca_id - 1;
			pdns_ca_C = pdns_ca_id + 1;
		}

		std::array<float ,3> pdns_ca_C2N_vec = getVector(xyzs, pdns_ca_N) - getVector(xyzs, pdns_ca_C);
		float pdns_ca_C2N_scholar = calc_Norm(pdns_ca_C2N_vec);
		std::array<float, 3> pdns_ca_vec = getVector(xyzs, pdns_ca_id);

		for (const std::array<int, 2>& phos_sug_id : pdns_phos_sug_ids) {
			std::array<float, 3> pdns_phos_vec = getVector(xyzs, phos_sug_id[0]);
			if (calcDistance(pdns_phos_vec, pdns_ca_vec) > (5.0 + parameter[0])) continue;

			std::array<float, 3> pdns_phos2ca_vec = pdns_ca_vec - pdns_phos_vec;
			float pdns_phos2ca_scholar = calc_Norm(pdns_phos2ca_vec);
			std::array<float, 3> pdns_phos2sug_vec = getVector(xyzs, phos_sug_id[1]) - pdns_phos_vec;
			float pdns_phos2sug_scholar = calc_Norm(pdns_phos2sug_vec);

			float cos_ca2phos2sug = calcDot(pdns_phos2ca_vec, pdns_phos2sug_vec) / (pdns_phos2ca_scholar * pdns_phos2sug_scholar);
			float cos_N2C2phos = calcDot(pdns_ca_C2N_vec, pdns_phos2ca_vec) / (pdns_ca_C2N_scholar * pdns_phos2ca_scholar);
			float rad_ca2phos2sug = acos(cos_ca2phos2sug);
			float rad_N2C2phos = acos(cos_N2C2phos);

			float d_rad_ca2phos2sug = rad_ca2phos2sug - parameter[1];
			float d_rad_N2C2phos = rad_N2C2phos - parameter[2];
			if (((-2 * cutoff_angle) <= d_rad_ca2phos2sug) &&
				(d_rad_ca2phos2sug <= 2 * cutoff_angle) &&
				((-2 * cutoff_angle) <= d_rad_N2C2phos) &&
				(d_rad_N2C2phos <= 2 * cutoff_angle)) {
				result.push_back({phos_sug_id[0], pdns_ca_id});
				break;
			}
		}
	}
	return result;
}


std::vector<std::vector<int>> cafemol::DCDAnalyzer::get_PDNSContactResidues(const std::vector<int>& pdns_pro_ids, const std::vector<int>& dna_ids, const float& cutoff) {

    read_num_frame_and_atom();
    std::vector<std::vector<int>> result(frame_num);
    for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
        std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
        std::vector<int> pdns_target_ids;
        for (const int& pdns_pro_id : pdns_pro_ids) {
            std::array<float, 3> pdns_pro_vec = getVector(xyz, pdns_pro_id);
            pdns_target_ids.push_back(getClosestResidue(xyz, pdns_pro_vec, dna_ids, cutoff));
        }
        result[iframe] = pdns_target_ids;
    }
    return result;
}


std::vector<std::vector<std::array<int, 2>>> cafemol::DCDAnalyzer::get_PDNSContactResidues(const std::vector<std::array<int, 2>>& pdns_phos_sug_ids, const std::vector<int>& pdns_pro_ids, const std::vector<std::array<float ,3>>& parameters) {

	read_num_frame_and_atom();
	std::vector<std::vector<std::array<int, 2>>> result(frame_num);
	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		result[iframe] = get_PDNSContactResidues(xyz, pdns_phos_sug_ids, pdns_pro_ids, parameters);
	}
	load_file();
	return result;
}


std::vector<bool> cafemol::DCDAnalyzer::is_NativeContacts(const std::array<std::vector<float>, 3>& xyzs, const std::vector<int>& contact_search_ids, const std::vector<int>& contact_target_ids, const float& cutoff) {

    std::vector<bool> result;

    for (const int& contact_search_id : contact_search_ids) {
        std::array<float, 3> contact_search_vec = getVector(xyzs, contact_search_id);
        std::array<float, 3> contact_target_vec;
        bool is_contacted = false;
        for (const int& contact_target_id : contact_target_ids) {
            contact_target_vec = getVector(xyzs, contact_target_id);
            if (is_NativeContact(contact_search_vec, contact_target_vec, cutoff)) {
                is_contacted = true;
                break;
            }
        }
        result.push_back(is_contacted);
    }
    return result;
}


std::vector<std::vector<bool>> cafemol::DCDAnalyzer::is_NativeContacts(const std::vector<int>& contact_search_ids, const std::vector<int>& contact_target_ids, const float& cutoff) {

    read_num_frame_and_atom();
    std::vector<std::vector<bool>> result(frame_num);

    for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
        std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
        result[iframe] = is_NativeContacts(xyz, contact_search_ids, contact_target_ids, cutoff);
    }
    load_file();
    return result;
}
