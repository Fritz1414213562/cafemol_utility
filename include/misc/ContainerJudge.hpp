#ifndef CONTAINER_JUDGE_HPP
#define CONTAINER_JUDGE_HPP
#include<algorithm>
#include<vector>


namespace cafemol {
	
namespace library {


template<typename Type>
bool is_contains(const std::vector<Type>& container, const Type& component) {
	bool result = false;
	for (const Type& container_component : container) {
		result = result || (container_component == component);
	}
	return result;
}


template<typename Container>
bool is_containercontains(const std::vector<Container>& containers, const Container& container) {
	bool result = false;
	for (const Container& comp_container : containers) {
		result = result ||
				 (comp_container.size() == container.size() &&
				 std::equal(comp_container.cbegin(), comp_container.cend(), container.cbegin()));
	}
	return result;
}


}
}

#endif /* CONTAINER_JUDGE_HPP */
