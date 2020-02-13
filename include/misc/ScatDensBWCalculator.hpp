#ifndef SCAT_DENSBW_CALCULATOR_HPP
#define SCAT_DENSBW_CALCULATOR_HPP
#include<UtilFunc.hpp>
#include<SetSimilarityCalculation.hpp>
#include<OtherFormat/ClusteringResultReader.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<vector>
#include<array>
#include<string>
#include<cmath>
#include<memory>


namespace cafemol {

namespace analysis {


enum Scat_DensBW_Dist_Frags {

	SCAT_DENSBW_DIST_SS,
	SCAT_DENSBW_DIST_JS,
	SCAT_DENSBW_DIST_DS,

};



class S_DbwCalculator {

public:
	
	S_DbwCalculator(const float& bin_w) : bin_width(bin_w) {}

	~S_DbwCalculator() = default;


	template<Scat_DensBW_Dist_Frags Dist_Frag>
	float run(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::vector<std::vector<std::size_t>>& cluster_indices_set) {
		return calc_S_Dbw<Dist_Frag>(all_data, cluster_indices_set);
	}


	template<Scat_DensBW_Dist_Frags Dist_Frag>
	void run(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::vector<std::vector<std::size_t>>& cluster_indices_set,
	std::vector<float>& out) {
		 out.push_back(calc_S_Dbw<Dist_Frag>(all_data, cluster_indices_set));
	}


	template<Scat_DensBW_Dist_Frags Dist_Frag>
	float run(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::string& clustering_result_filename) {

		std::unique_ptr<cafemol::ClusteringResultReader> clustering_result_reader = std::make_unique<cafemol::ClusteringResultReader>(clustering_result_filename);

		const std::vector<std::vector<std::size_t>>& cluster_indices_set = clustering_result_reader->read_ClassifiedIDs();

		sout("Read the file, " + clustering_result_filename);
		return calc_S_Dbw<Dist_Frag>(all_data, cluster_indices_set);
	}


	template<Scat_DensBW_Dist_Frags Dist_Frag>
	void run(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::string& clustering_result_filename,
	std::vector<float>& out) {

		std::unique_ptr<cafemol::ClusteringResultReader> clustering_result_reader = std::make_unique<cafemol::ClusteringResultReader>(clustering_result_filename);

		const std::vector<std::vector<std::size_t>>& cluster_indices_set = clustering_result_reader->read_ClassifiedIDs();

		sout("Read the file, " + clustering_result_filename);
		out.push_back(calc_S_Dbw<Dist_Frag>(all_data, cluster_indices_set));
	}



private:


// member variants

	const float bin_width = 1.0;
	const std::size_t BLOCK_SIZE = 90;


	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();


// member functions



	template<Scat_DensBW_Dist_Frags Dist_Frag>
	const float calc_S_Dbw(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::vector<std::vector<std::size_t>>& cluster_indices_set) {
		
		// S_Dbw(k) = Scat(k) + Dens_bw(k)

		sout.output_HyphenBlock("", BLOCK_SIZE);
		sout("Calculation starts.");

		// calculate centroids
		const std::size_t& cluster_number = cluster_indices_set.size();
		if (cluster_number <= 1) {
			eout(
			"In 'calc_S_Dbw'",
			"The cluster number is less than 1. It must be more than 2");
		}
		const float& f_cluster_number = static_cast<float>(cluster_number);

		std::vector<std::vector<std::array<float, 3>>> centroids;
		for (const std::vector<std::size_t>& cluster_indices : cluster_indices_set) {
			centroids.push_back(calc_ClusterCentroid(all_data, cluster_indices));
		}



		// calculate Scat(k)

		sout(
		"",
		"calculating Scatter Extent 'Scat(k)'");
		
		const std::vector<float>& variances = calc_VariancesOfClusters<Dist_Frag>(all_data, cluster_indices_set, centroids);
		const float& sse_of_all = calc_SumSquaredErrorOfAllData<Dist_Frag>(all_data);
		float scatter_extent = 0.0;

		for (std::size_t i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			scatter_extent += (
			variances[i_cluster] / (sse_of_all / static_cast<float>(cluster_indices_set[i_cluster].size())));
		}

		scatter_extent /= f_cluster_number;

	//	sout(scatter_extent);


		sout(".... Done");

		// calculate Dens_bw(k)

		sout(
		"",
		"calculating Density Extent 'Dens_bw(k)'");


		float sum_of_variances = 0.0;

		for (const float& variance_of_cluster : variances) {
			sum_of_variances += variance_of_cluster;
		}
		const float& standard_deviation = std::sqrt(sum_of_variances) / f_cluster_number;
		sout(standard_deviation);

		float density_extent = 0.0;

		for (std::size_t i_cluster = 0; i_cluster < cluster_number - 1; ++i_cluster) {

			// count data number near cluster i

			const std::vector<std::array<float, 3>>& centroid_i = centroids[i_cluster];
			const std::vector<std::size_t>& cluster_indices_i = cluster_indices_set[i_cluster];
			const float& data_num_near_i = static_cast<float>(count_DataNearCentroid<Dist_Frag>(all_data, cluster_indices_i, centroid_i, standard_deviation));


			for (std::size_t j_cluster = i_cluster + 1; j_cluster < cluster_number; ++j_cluster) {
				// count data number near cluster j
				const std::vector<std::array<float, 3>>& centroid_j = centroids[j_cluster];
				const std::vector<std::size_t>& cluster_indices_j = cluster_indices_set[j_cluster];
				const float& data_num_near_j = static_cast<float>(count_DataNearCentroid<Dist_Frag>(all_data, cluster_indices_j, centroid_j, standard_deviation));

				// count data number near the mid point between the centroids of cluster i and j

			//	const std::vector<std::array<float, 3>>& centroid_of_ij_union = calc_ClusterCentroid(all_data, combined_ij_indices);

				const std::vector<std::size_t>& combined_ij_indices = combine_ClusterIndices(cluster_indices_i, cluster_indices_j);
			//	const float& data_num_near_ij_midium = static_cast<float>(count_DataNearCentroid(all_data, combined_ij_indices, centroid_of_ij_union, standard_deviation));

				std::vector<std::array<float, 3>> centroid_of_ij_union = centroid_i;
				cafemol::library::overlie_Data(centroid_of_ij_union, centroid_j, bin_width);

				for (std::size_t data_index = 0; data_index < centroid_of_ij_union.size(); ++data_index) {
					centroid_of_ij_union[data_index][2] /= 2.0;
				}

				const float& data_num_near_ij_midium = static_cast<float>(count_DataNearCentroid<Dist_Frag>(all_data, combined_ij_indices, centroid_of_ij_union, standard_deviation));

				// calculate density

				const float& max_data_num_bw_ij = std::max(data_num_near_i, data_num_near_j);

				density_extent += data_num_near_ij_midium / max_data_num_bw_ij;

			}
		}

		density_extent *= (2 / (f_cluster_number * (f_cluster_number - 1.)));


		sout(density_extent);
		sout(".... Done");

		sout(
		"Calculation ends",
		"");


		return scatter_extent + density_extent;

	}







// sse calculation methods



	template<Scat_DensBW_Dist_Frags Dist_Frag>
	const std::vector<float> calc_VariancesOfClusters(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::vector<std::vector<std::size_t>>& cluster_indices_set,
	const std::vector<std::vector<std::array<float, 3>>>& centroids) {

		if (cluster_indices_set.size() != centroids.size()) {
			eout(
			"In 'calc_VariancesOfClusters'", 
			"The size of cluster indices set must be consistent with that of centoids");
		}


		const std::size_t& cluster_number = cluster_indices_set.size();
		std::vector<float> result(cluster_number, 0);

		for (std::size_t cluster_id = 0; cluster_id < cluster_number; ++cluster_id) {
			const std::vector<std::size_t>& cluster_indices = cluster_indices_set[cluster_id];
			const std::vector<std::array<float, 3>>& centroid = centroids[cluster_id];

			result[cluster_id] = calc_Variance<Dist_Frag>(all_data, cluster_indices, centroid);
		}

		return result;
	}





	template<Scat_DensBW_Dist_Frags Dist_Frag>
	const float calc_Variance(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::vector<std::size_t>& cluster_indices,
	const std::vector<std::array<float, 3>>& centroid) {

		const float& cluster_size = static_cast<float>(cluster_indices.size());

		float sum_of_squared_error = 0.0;

		for (const std::size_t& cluster_index : cluster_indices) {
			const float& similarity_distance = calc_SimilarityDistance<Dist_Frag>(all_data[cluster_index], centroid);
			sum_of_squared_error += (similarity_distance * similarity_distance);
		}

		const float& result = sum_of_squared_error / cluster_size;

		return result;
	}






	template<Scat_DensBW_Dist_Frags Dist_Frag>
	const float calc_SumSquaredErrorOfAllData(const std::vector<std::vector<std::array<float, 3>>>& all_data) {
		
		float result = 1.0;

		const std::vector<std::array<float, 3>>& centroid = calc_CentroidFromAll(all_data);

		float sum_of_squared_error = 0.0;

		for (const std::vector<std::array<float, 3>>& frame : all_data) {
			const float& similarity_distance = calc_SimilarityDistance<Dist_Frag>(centroid, frame);
			sum_of_squared_error += (similarity_distance * similarity_distance);
		}
		result = sum_of_squared_error;

		return result;

	}




// density calculation methods



	template<Scat_DensBW_Dist_Frags Dist_Frag>
	const std::size_t count_DataNearCentroid(
	const std::vector<std::vector<std::array<float, 3>>>& all_data,
	const std::vector<std::size_t>& cluster_indices,
	const std::vector<std::array<float, 3>>& centroid,
	const float& std_dev) {
		std::size_t result = 0;

		for (const std::size_t& cluster_index : cluster_indices) {
			const std::vector<std::array<float, 3>>& frame = all_data[cluster_index];
			result += static_cast<std::size_t>(is_Near2Centroid<Dist_Frag>(frame, centroid, std_dev));
		}

		return result;
	}





	template<Scat_DensBW_Dist_Frags Dist_Frag>
	const bool is_Near2Centroid(const std::vector<std::array<float, 3>>& frame, const std::vector<std::array<float, 3>>& centroid, const float& std_dev) {
		return (calc_SimilarityDistance<Dist_Frag>(frame, centroid) <= std_dev);
	}




// common use in this class



	const std::vector<std::array<float, 3>> calc_ClusterCentroid(const std::vector<std::vector<std::array<float, 3>>>& all_data, const std::vector<std::size_t>& cluster_indices) {

		std::vector<std::array<float, 3>> result = all_data[cluster_indices[0]];
		const std::size_t& cluster_size = cluster_indices.size();

		for (std::size_t index = 1; index < cluster_size; ++index) {
			const std::size_t& cluster_index = cluster_indices[index];
			cafemol::library::overlie_Data(result, all_data[cluster_index], bin_width);
		}

		const std::size_t& data_size = result.size();
		for (std::size_t data_index = 0; data_index < data_size; ++data_index) {
			result[data_index][2] /= static_cast<float>(cluster_size);
		}

		return result;
	}



	const std::vector<std::array<float, 3>> calc_CentroidFromAll(const std::vector<std::vector<std::array<float, 3>>>& all_data) {

		std::vector<std::array<float, 3>> result = all_data[0];
		const std::size_t& frame_size = all_data.size();

		for (std::size_t index = 1; index < frame_size; ++index) {
			cafemol::library::overlie_Data(result, all_data[index], bin_width);
		}

		const std::size_t& data_size = result.size();
		for (std::size_t data_index = 0; data_index < data_size; ++data_index) {
			result[data_index][2] /= static_cast<float>(frame_size);
		}

		return result;
	}



	const std::vector<std::size_t> combine_ClusterIndices(const std::vector<std::size_t>& lhs, const std::vector<std::size_t>& rhs) {
		std::vector<std::size_t> result;

		for (const std::size_t& lhs_index : lhs) result.push_back(lhs_index);
		for (const std::size_t& rhs_index : rhs) result.push_back(rhs_index);

		return result;
	}


	// instanciated bottom this file
	template<Scat_DensBW_Dist_Frags Dist_Frag>
	const float calc_SimilarityDistance(const std::vector<std::array<float, 3>>& lhs, const std::vector<std::array<float, 3>>& rhs);


};




/*


namespace cafemol::analysis


*/


// similarity calculation methods


template<>
inline const float
cafemol::analysis::S_DbwCalculator::calc_SimilarityDistance<
cafemol::analysis::SCAT_DENSBW_DIST_SS>(
const std::vector<std::array<float, 3>>& lhs, const std::vector<std::array<float, 3>>& rhs) {
	return cafemol::library::calc_SimpsonSimilarity(lhs, rhs, bin_width);
}



template<>
inline const float
cafemol::analysis::S_DbwCalculator::calc_SimilarityDistance<
cafemol::analysis::SCAT_DENSBW_DIST_JS>(
const std::vector<std::array<float, 3>>& lhs, const std::vector<std::array<float, 3>>& rhs) {
	return cafemol::library::calc_JaccardSimilarity(lhs, rhs, bin_width);
}



template<>
inline const float
cafemol::analysis::S_DbwCalculator::calc_SimilarityDistance<
cafemol::analysis::SCAT_DENSBW_DIST_DS>(
const std::vector<std::array<float, 3>>& lhs, const std::vector<std::array<float, 3>>& rhs) {
	return cafemol::library::calc_DiceSimilarity(lhs, rhs, bin_width);
}



}
} // terminate cafemol::analysis scope


#endif /*SCAT_DENSBW_CALCULATOR_HPP*/
