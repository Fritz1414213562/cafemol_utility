#ifndef CONTACT_STATE_CLUSTERING_HPP
#define CONTACT_STATE_CLUSTERING_HPP
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<opencv2/core.hpp>
//#include<opencv2/imgproc.hpp>
#include<SetSimilarityCalculation.hpp>
#include<SetDistanceCalculation.hpp>
#include<UtilFunc.hpp>
#include<random>
#include<iostream>
#include<string>
#include<vector>
#include<array>
#include<tuple>
#include<cmath>


//namespace cafemol {
//
//namespace analysis {

namespace cafemol {
	using ClusteringResults = std::tuple<std::vector<std::vector<std::size_t>>,
										 std::vector<float>,
										 std::vector<std::vector<float>>>;
}

namespace cafemol::analysis {

enum KMeans_Frags {

// distance
	// EMD
	KMF_DIST_EMD,

	// Szymkiewicz-Simpson Similarity
	KMF_DIST_SS,

	// Dice Similarity
	KMF_DIST_DS,

	// Jaccard Similarity
	KMF_DIST_JS,

// initialization
	KMF_INI_DEFAULT,
	KMF_INI_PLUSPLUS,

};

}

namespace cafemol::analysis {

class ContactStateClustering {

public:


	// constructor

	ContactStateClustering() {
		mt_engine.seed(random_seed);
	}

	ContactStateClustering(const std::size_t& cluster_num)
		: n_cluster(cluster_num) {
			mt_engine.seed(random_seed);
	}

	ContactStateClustering(const std::size_t& cluster_num, const std::size_t& iteration_num)
		: n_cluster(cluster_num), MAX_ITERATATION(iteration_num) {
			mt_engine.seed(random_seed);
	}

	ContactStateClustering(const std::size_t& cluster_num, const std::size_t& iteration_num, const float& cutoff_siml)
		: n_cluster(cluster_num), MAX_ITERATATION(iteration_num), cutoff_similarity(cutoff_siml) {
			mt_engine.seed(random_seed);
	}

	ContactStateClustering(const std::size_t& cluster_num, const std::size_t& iteration_num, const float& cutoff_siml, const std::size_t& n_seed)
		: n_cluster(cluster_num), MAX_ITERATATION(iteration_num), cutoff_similarity(cutoff_siml), random_seed(n_seed) {
			mt_engine.seed(random_seed);
	}

	ContactStateClustering(const std::size_t& cluster_num, const std::size_t& iteration_num, const float& cutoff_siml, const float& bin_w)
		: n_cluster(cluster_num), MAX_ITERATATION(iteration_num), cutoff_similarity(cutoff_siml), bin_width(bin_w) {
			mt_engine.seed(random_seed);
	}

	ContactStateClustering(const std::size_t& cluster_num, const std::size_t& iteration_num, const float& cutoff_siml, const std::size_t& n_seed, const float& bin_w)
		: n_cluster(cluster_num), MAX_ITERATATION(iteration_num), cutoff_similarity(cutoff_siml), random_seed(n_seed), bin_width(bin_w) {
			mt_engine.seed(random_seed);
	}



	// deconstructor
	~ContactStateClustering() = default;


	// run


	template<typename Size3_Vector, KMeans_Frags Dist_Frag, KMeans_Frags Ini_Frag>
	const std::tuple<std::vector<std::vector<std::size_t>>,
					 std::vector<float>,
					 std::vector<std::vector<float>>>
	run(const std::vector<std::vector<Size3_Vector>>& all_data) {


		sout("make initial cluster nuclears");

		std::vector<std::vector<std::size_t>> clusters_indices = init_ClusterNuclear<Size3_Vector, Dist_Frag, Ini_Frag>(all_data);
		sout("calculate initial centroids", "");
		std::vector<std::vector<Size3_Vector>> centroids = calc_Centroids(all_data,
																		  clusters_indices);
		std::tuple<std::vector<std::vector<std::size_t>>,
				   std::vector<float>,
				   std::vector<std::vector<float>>> result;
		std::get<1>(result).resize(n_cluster);
		std::get<2>(result).resize(n_cluster);
		for (std::size_t i_cluster = 0; i_cluster < n_cluster; ++i_cluster) {
			std::get<2>(result)[i_cluster].resize(n_cluster);
		}
		bool is_to_succeed = false;


		sout("Iteration starts.");
		std::size_t iterated_num = 1;


		for (std::size_t i_update = 0; i_update < MAX_ITERATATION; ++i_update) {
			const std::tuple<std::vector<std::vector<std::size_t>>, std::vector<float>>& clusters_indices_and_dist = classify_AllDataIntoClusters<Size3_Vector, Dist_Frag>(all_data, centroids);
			const std::vector<std::vector<Size3_Vector>>&
			updated_centroids = calc_Centroids(all_data, std::get<0>(clusters_indices_and_dist));
			const std::vector<float>& centroids_dist = measure_CentroidDistance<Size3_Vector, Dist_Frag>(centroids, updated_centroids);
		//	if (is_similar20vec<Size3_Vector>(centroids_dist)) {
			if (cafemol::library::is_similar20vec(centroids_dist, cutoff_similarity)) {
				std::get<0>(result) = std::get<0>(clusters_indices_and_dist);
				std::get<1>(result) = std::get<1>(clusters_indices_and_dist);
				for (std::size_t i_cluster = 0; i_cluster < n_cluster - 1; ++i_cluster) {
					for (std::size_t j_cluster = i_cluster; j_cluster < n_cluster; ++j_cluster) {
						const float& diff_centroids = calc_DataDistance<Size3_Vector, Dist_Frag>(updated_centroids[i_cluster], updated_centroids[j_cluster]);
						if (i_cluster != j_cluster) {
							std::get<2>(result)[i_cluster][j_cluster] = diff_centroids;
							std::get<2>(result)[j_cluster][i_cluster] = diff_centroids;
						}
						else std::get<2>(result)[i_cluster][j_cluster] = diff_centroids;
					}
					if (i_cluster == n_cluster - 2) {
						const float& diff_centroids = calc_DataDistance<Size3_Vector, Dist_Frag>(updated_centroids[n_cluster - 1], updated_centroids[n_cluster - 1]);
						std::get<2>(result)[n_cluster - 1][n_cluster - 1] = diff_centroids;
					}
				}
				is_to_succeed = true;
				break;
			}
			else {
				centroids = updated_centroids;
			}
			++iterated_num;
		}

		if (!is_to_succeed) {
			eout("Iteration is finished, but clustering doesn't succeed. You should cycle more.");
		}
		else {
			sout("Iterated number is " + std::to_string(iterated_num));
		}

		return result;
	}





private:

	// default
	const std::size_t n_cluster = 3;
	const std::size_t MAX_ITERATATION = 100;
	const float cutoff_similarity = 1.0;
	const std::size_t random_seed = 10000000;
	const float bin_width = 1.0;
	std::mt19937_64 mt_engine;


	// error output
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	// standard output
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();



	// check whether the cluster size is 0 before using this method
	template<typename Size3_Vector>
	const std::vector<Size3_Vector> calc_ClusterCentroid(const std::vector<std::vector<Size3_Vector>>& all_data, const std::vector<std::size_t>& indices_in_a_cluster) {

		std::vector<Size3_Vector> result = all_data[indices_in_a_cluster[0]];
		const std::size_t& cluster_size = indices_in_a_cluster.size();

		for (std::size_t idx = 1; idx < cluster_size; ++idx) {
			const std::size_t& index_in_a_cluster = indices_in_a_cluster[idx];
			cafemol::library::overlie_Data(result, all_data[index_in_a_cluster], bin_width);
		}

		const std::size_t& data_size = result.size();
		for (std::size_t data_index = 0; data_index < data_size; ++data_index) {
			result[data_index][2] /= static_cast<float>(cluster_size);
		}

		return result;
	}



	template<typename Size3_Vector>
	const std::vector<std::vector<Size3_Vector>>
	calc_Centroids(const std::vector<std::vector<Size3_Vector>>& all_data,
	const std::vector<std::vector<std::size_t>>& cluster_indices) {

		std::vector<std::vector<Size3_Vector>> result(n_cluster);

		for (std::size_t i_cluster = 0; i_cluster < n_cluster; ++i_cluster) {
			if (cluster_indices[i_cluster].size() < 1) continue;
			result[i_cluster] = calc_ClusterCentroid(all_data, cluster_indices[i_cluster]);
		}
		return result;
	}


	template<typename Size3_Vector, KMeans_Frags Dist_Frag>
	const std::tuple<std::vector<std::vector<std::size_t>>, std::vector<float>>
	classify_AllDataIntoClusters(const std::vector<std::vector<Size3_Vector>>& all_data, 
	const std::vector<std::vector<Size3_Vector>>& centroids) {

		std::tuple<std::vector<std::vector<std::size_t>>, std::vector<float>> result;
		std::get<0>(result).resize(n_cluster);
		std::get<1>(result).resize(n_cluster);
	//	std::vector<float> intracluster_dist(n_cluster, 0.0);
		std::vector<float> intracluster_ses(n_cluster, 0.0);

		for (std::size_t data_index = 0; data_index < all_data.size(); ++data_index) {
			const std::vector<Size3_Vector>& frame = all_data[data_index];
			const std::tuple<std::size_t, float>& cluster_id_and_dist = calc_ClusterIDAndDistFromCentroid<Size3_Vector, Dist_Frag>(frame, centroids);
			std::get<0>(result)[std::get<0>(cluster_id_and_dist)].push_back(data_index);
//			if (intracluster_dist[std::get<0>(cluster_id_and_dist)] < std::get<1>(cluster_id_and_dist)) {
//				intracluster_dist[std::get<0>(cluster_id_and_dist)] = std::get<1>(cluster_id_and_dist);
//			}
			intracluster_ses[std::get<0>(cluster_id_and_dist)] += std::get<1>(cluster_id_and_dist);
		}

//		std::get<1>(result) = intracluster_dist;
		std::get<1>(result) = intracluster_ses;

		return result;
	}




	template<typename Size3_Vector, KMeans_Frags Dist_Frag>
	const std::vector<float>
	measure_CentroidDistance(const std::vector<std::vector<Size3_Vector>>& new_centroids, const std::vector<std::vector<Size3_Vector>>& old_centroids) {
	
		std::vector<float> result(n_cluster);
	
		for (std::size_t i_cluster = 0; i_cluster < n_cluster; ++i_cluster) {
			const std::vector<Size3_Vector>& new_centroid = new_centroids[i_cluster];
			const std::vector<Size3_Vector>& old_centroid = old_centroids[i_cluster];
			const float& distance_between_new_old = calc_DataDistance<Size3_Vector, Dist_Frag>(new_centroid, old_centroid);
			result[i_cluster] = distance_between_new_old;
		}
		return result;
	}



	template<typename Size3_Vector, KMeans_Frags Dist_Frag, KMeans_Frags Ini_Frag>
	const std::vector<std::vector<std::size_t>> init_ClusterNuclear(const std::vector<std::vector<Size3_Vector>>& all_data);


	
	template<typename Size3_Vector, KMeans_Frags Ini_Frag>
	const std::vector<std::vector<std::size_t>> init_ClusterNuclear(const std::vector<std::vector<Size3_Vector>>& all_data);


	
	template<typename Size3_Vector, KMeans_Frags Dist_Frag>
	const float calc_DataDistance(const std::vector<Size3_Vector>& lhs, const std::vector<Size3_Vector>& rhs);



	template<typename Size3_Vector, KMeans_Frags Dist_Frag>
	const std::tuple<std::size_t, float>
	calc_ClusterIDAndDistFromCentroid(const std::vector<Size3_Vector>& frame, const std::vector<std::vector<Size3_Vector>>& centroids);




// --------------------------------------------------------------------------------------
// Initial cluster-nuclear sampling method


	// totaly random
	const std::vector<std::vector<std::size_t>> choose_TotalyRandomIndices(const std::size_t& all_data_size) {
	
		std::mt19937_64 mt_engine(random_seed);
		std::uniform_int_distribution<> distribution(0, all_data_size - 1);
		std::vector<std::vector<std::size_t>> result(n_cluster);
	
		for (std::size_t i_cluster = 0; i_cluster < n_cluster; ++i_cluster) {
			result[i_cluster].push_back(distribution(mt_engine));
		}
	
		return result;
	}



	// using ++method (just choosing one cluster nuclear randomly, but others are chosen by the distance)

	template<typename Size3_Vector, KMeans_Frags Dist_Frag>
	const std::vector<std::vector<std::size_t>> choose_DistantIndices(const std::vector<std::vector<Size3_Vector>>& all_data) {

		std::mt19937_64 mt_engine_for_uniform_dist(random_seed);
		std::mt19937_64 mt_engine_for_piecewise_dist(random_seed + 1);
		std::uniform_int_distribution<> uni_int_distribution(0, all_data.size() - 1);
		std::vector<std::vector<std::size_t>> result(n_cluster);

		result[0].push_back(uni_int_distribution(mt_engine_for_uniform_dist));

		for (std::size_t i_cluster = 1; i_cluster < n_cluster; ++i_cluster) {
			double total_dist = 0.0;
			std::vector<double> prob_intervals(all_data.size() + 1, 0);
			std::vector<double> prob_densities(all_data.size());
			std::size_t datum_index = 0;
			for (const std::vector<Size3_Vector>& datum : all_data) {
				double dist = 0.0;
				for (std::size_t j_cluster = 0; j_cluster < i_cluster; ++j_cluster) {
					const float& data_dist_between_j_datum = calc_DataDistance<Size3_Vector, Dist_Frag>(datum, all_data[result[j_cluster][0]]);
					dist += static_cast<double>(data_dist_between_j_datum);
					total_dist += static_cast<double>(data_dist_between_j_datum);
				}
				prob_intervals[datum_index + 1] = static_cast<double>(datum_index + 1);
				prob_densities[datum_index] = dist;
				++datum_index;
			}

			for (std::size_t prob_dens_idx = 0; prob_dens_idx < prob_densities.size(); ++prob_dens_idx) {
				prob_densities[prob_dens_idx] /= total_dist;
			}

			std::piecewise_constant_distribution<> pw_const_dist(
				prob_intervals.begin(),
				prob_intervals.end(),
				prob_densities.begin()
			);
			
			const std::size_t& nuclear_index = std::floor(pw_const_dist(mt_engine_for_piecewise_dist));
			result[i_cluster].push_back(nuclear_index);
		}

		return result;
	}






// --------------------------------------------------------------------------------------
// Distance Calculation Methods





	// Szymkiewicz-Simpson coefficient

	template<typename Size3_Vector>
	const float calculate_SimpsonSimilarity(const std::vector<Size3_Vector>& lhs, const std::vector<Size3_Vector>& rhs) {
		return cafemol::library::calc_SimpsonSimilarity(lhs, rhs, bin_width);
	}

	// Dice Index

	template<typename Size3_Vector>
	const float calculate_DiceSimilarity(const std::vector<Size3_Vector>& lhs, const std::vector<Size3_Vector>& rhs) {
		return cafemol::library::calc_DiceSimilarity(lhs, rhs, bin_width);
	}

	// Jaccard Index

	template<typename Size3_Vector>
	const float calculate_JaccardSimilarity(const std::vector<Size3_Vector>& lhs, const std::vector<Size3_Vector>& rhs) {
		return cafemol::library::calc_JaccardSimilarity(lhs, rhs, bin_width);
	}




// -------------------------------------------------------------------------------------------
// Data Classification Methods


	// using the real distance between frame and centroids
	template<typename Size3_Vector, KMeans_Frags Dist_Frag>
	const std::tuple<std::size_t, float> 
	calc_ClosestCentroidIDAndDist(const std::vector<Size3_Vector>& frame, const std::vector<std::vector<Size3_Vector>>& centroids) {
	
		const std::vector<Size3_Vector>& centroid_index0 = centroids[0];
		float min_data_distance = calc_DataDistance<Size3_Vector, Dist_Frag>(frame, centroid_index0);
		std::tuple<std::size_t, float> result = std::make_tuple(0, min_data_distance);
		for (std::size_t i_cluster = 1; i_cluster < n_cluster; ++i_cluster) {
			float data_distance = calc_DataDistance<Size3_Vector, Dist_Frag>(frame, centroids[i_cluster]);
			if (data_distance < min_data_distance) {
				min_data_distance = data_distance;
				std::get<0>(result) = i_cluster;
				std::get<1>(result) = min_data_distance;
			}
		}
		return result;
	}


	// using the similarity between frame and centroids

	template<typename Size3_Vector, KMeans_Frags Dist_Frag>
	const std::tuple<std::size_t, float>
	calc_MostSimilarCentroidIDAndDist(const std::vector<Size3_Vector>& frame, const std::vector<std::vector<Size3_Vector>>& centroids) {

		const std::vector<Size3_Vector>& centroid_index0 = centroids[0];
		float min_data_similarity = calc_DataDistance<Size3_Vector, Dist_Frag>(frame, centroid_index0);
		std::tuple<std::size_t, float> result = std::make_tuple(0, min_data_similarity);
		for (std::size_t i_cluster = 1; i_cluster < n_cluster; ++i_cluster) {
			float data_similarity = calc_DataDistance<Size3_Vector, Dist_Frag>(frame, centroids[i_cluster]);
			if (data_similarity < min_data_similarity) {
				min_data_similarity = data_similarity;
				std::get<0>(result) = i_cluster;
				std::get<1>(result) = min_data_similarity;
			}
		}
		return result;
	}



};

// -------------------------------------------------------------------------------------------
/* Global Scope



*/





// ------------------------------------------------------------------------------------------
// template instantiation



// initialization of the nuclear of clusters

// totaly random

template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<cv::Vec3f, cafemol::analysis::KMF_INI_DEFAULT>(const std::vector<std::vector<cv::Vec3f>>& all_data) {return choose_TotalyRandomIndices(all_data.size());}



template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<std::array<float, 3>, cafemol::analysis::KMF_INI_DEFAULT>(const std::vector<std::vector<std::array<float, 3>>>& all_data) {return choose_TotalyRandomIndices(all_data.size());}



// efficient sampling

	// EMD

template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<cv::Vec3f, cafemol::analysis::KMF_DIST_EMD, cafemol::analysis::KMF_INI_PLUSPLUS>(const std::vector<std::vector<cv::Vec3f>>& all_data) {return choose_DistantIndices<cv::Vec3f, cafemol::analysis::KMF_DIST_EMD>(all_data);}


	// SS

template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<cv::Vec3f, cafemol::analysis::KMF_DIST_SS, cafemol::analysis::KMF_INI_PLUSPLUS>(const std::vector<std::vector<cv::Vec3f>>& all_data) {return choose_DistantIndices<cv::Vec3f, cafemol::analysis::KMF_DIST_SS>(all_data);}



template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<std::array<float, 3>, cafemol::analysis::KMF_DIST_SS, cafemol::analysis::KMF_INI_PLUSPLUS>(const std::vector<std::vector<std::array<float, 3>>>& all_data) {return choose_DistantIndices<std::array<float, 3>, cafemol::analysis::KMF_DIST_SS>(all_data);}




	// JS

template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<cv::Vec3f, cafemol::analysis::KMF_DIST_JS, cafemol::analysis::KMF_INI_PLUSPLUS>(const std::vector<std::vector<cv::Vec3f>>& all_data) {return choose_DistantIndices<cv::Vec3f, cafemol::analysis::KMF_DIST_JS>(all_data);}



template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<std::array<float, 3>, cafemol::analysis::KMF_DIST_JS, cafemol::analysis::KMF_INI_PLUSPLUS>(const std::vector<std::vector<std::array<float, 3>>>& all_data) {return choose_DistantIndices<std::array<float, 3>, cafemol::analysis::KMF_DIST_JS>(all_data);}





	// DS

template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<cv::Vec3f, cafemol::analysis::KMF_DIST_DS, cafemol::analysis::KMF_INI_PLUSPLUS>(const std::vector<std::vector<cv::Vec3f>>& all_data) {return choose_DistantIndices<cv::Vec3f, cafemol::analysis::KMF_DIST_DS>(all_data);}



template<>
inline const std::vector<std::vector<std::size_t>> ContactStateClustering::init_ClusterNuclear<std::array<float, 3>, cafemol::analysis::KMF_DIST_DS, cafemol::analysis::KMF_INI_PLUSPLUS>(const std::vector<std::vector<std::array<float, 3>>>& all_data) {return choose_DistantIndices<std::array<float, 3>, cafemol::analysis::KMF_DIST_DS>(all_data);}






// measuring the distance among data


	// using EMD
template<>
inline const float
ContactStateClustering::calc_DataDistance<cv::Vec3f, cafemol::analysis::KMF_DIST_EMD>(
const std::vector<cv::Vec3f>& lhs, const std::vector<cv::Vec3f>& rhs) {
	return cafemol::library::calc_EarthMoversDistance(lhs, rhs);
}


	// using SS
template<>
inline const float ContactStateClustering::calc_DataDistance<cv::Vec3f, cafemol::analysis::KMF_DIST_SS>(const std::vector<cv::Vec3f>& lhs, const std::vector<cv::Vec3f>& rhs) {return calculate_SimpsonSimilarity(lhs, rhs);}



template<>
inline const float ContactStateClustering::calc_DataDistance<std::array<float, 3>, cafemol::analysis::KMF_DIST_SS>(const std::vector<std::array<float, 3>>& lhs, const std::vector<std::array<float, 3>>& rhs) {return calculate_SimpsonSimilarity(lhs, rhs);}



	// using JS
template<>
inline const float ContactStateClustering::calc_DataDistance<cv::Vec3f, cafemol::analysis::KMF_DIST_JS>(const std::vector<cv::Vec3f>& lhs, const std::vector<cv::Vec3f>& rhs) {return calculate_JaccardSimilarity(lhs, rhs);}



template<>
inline const float ContactStateClustering::calc_DataDistance<std::array<float, 3>, cafemol::analysis::KMF_DIST_JS>(const std::vector<std::array<float, 3>>& lhs, const std::vector<std::array<float, 3>>& rhs) {return calculate_JaccardSimilarity(lhs, rhs);}





	// using DS
template<>
inline const float ContactStateClustering::calc_DataDistance<cv::Vec3f, cafemol::analysis::KMF_DIST_DS>(const std::vector<cv::Vec3f>& lhs, const std::vector<cv::Vec3f>& rhs) {return calculate_DiceSimilarity(lhs, rhs);}



template<>
inline const float ContactStateClustering::calc_DataDistance<std::array<float, 3>, cafemol::analysis::KMF_DIST_DS>(const std::vector<std::array<float, 3>>& lhs, const std::vector<std::array<float, 3>>& rhs) {return calculate_DiceSimilarity(lhs, rhs);}





// searching the closest centroid
	// using real distance
template<>
inline const std::tuple<std::size_t, float> ContactStateClustering::calc_ClusterIDAndDistFromCentroid<cv::Vec3f, cafemol::analysis::KMF_DIST_EMD>(const std::vector<cv::Vec3f>& frame, const std::vector<std::vector<cv::Vec3f>>& centroids) {
	return calc_ClosestCentroidIDAndDist<cv::Vec3f, cafemol::analysis::KMF_DIST_EMD>(frame, centroids);
}


	// using simpson similarity
template<>
inline const std::tuple<std::size_t, float> ContactStateClustering::calc_ClusterIDAndDistFromCentroid<std::array<float, 3>, cafemol::analysis::KMF_DIST_SS>(const std::vector<std::array<float, 3>>& frame, const std::vector<std::vector<std::array<float, 3>>>& centroids) {return calc_MostSimilarCentroidIDAndDist<std::array<float, 3>, cafemol::analysis::KMF_DIST_SS>(frame, centroids);}


template<>
inline const std::tuple<std::size_t, float> ContactStateClustering::calc_ClusterIDAndDistFromCentroid<cv::Vec3f, cafemol::analysis::KMF_DIST_SS>(const std::vector<cv::Vec3f>& frame, const std::vector<std::vector<cv::Vec3f>>& centroids) {return calc_MostSimilarCentroidIDAndDist<cv::Vec3f, cafemol::analysis::KMF_DIST_SS>(frame, centroids);}



	// using jaccard similarity

template<>
inline const std::tuple<std::size_t, float> ContactStateClustering::calc_ClusterIDAndDistFromCentroid<std::array<float, 3>, cafemol::analysis::KMF_DIST_JS>(const std::vector<std::array<float, 3>>& frame, const std::vector<std::vector<std::array<float, 3>>>& centroids) {return calc_MostSimilarCentroidIDAndDist<std::array<float, 3>, cafemol::analysis::KMF_DIST_JS>(frame, centroids);}


template<>
inline const std::tuple<std::size_t, float> ContactStateClustering::calc_ClusterIDAndDistFromCentroid<cv::Vec3f, cafemol::analysis::KMF_DIST_JS>(const std::vector<cv::Vec3f>& frame, const std::vector<std::vector<cv::Vec3f>>& centroids) {return calc_MostSimilarCentroidIDAndDist<cv::Vec3f, cafemol::analysis::KMF_DIST_JS>(frame, centroids);}





	// using dice similarity

template<>
inline const std::tuple<std::size_t, float> ContactStateClustering::calc_ClusterIDAndDistFromCentroid<std::array<float, 3>, cafemol::analysis::KMF_DIST_DS>(const std::vector<std::array<float, 3>>& frame, const std::vector<std::vector<std::array<float, 3>>>& centroids) {return calc_MostSimilarCentroidIDAndDist<std::array<float, 3>, cafemol::analysis::KMF_DIST_DS>(frame, centroids);}


template<>
inline const std::tuple<std::size_t, float> ContactStateClustering::calc_ClusterIDAndDistFromCentroid<cv::Vec3f, cafemol::analysis::KMF_DIST_DS>(const std::vector<cv::Vec3f>& frame, const std::vector<std::vector<cv::Vec3f>>& centroids) {return calc_MostSimilarCentroidIDAndDist<cv::Vec3f, cafemol::analysis::KMF_DIST_DS>(frame, centroids);}







// terminate the namespace cafemol::analysis
}
//}
//}

#endif /* CONTACT_STATE_CLUSTERING_HPP */
