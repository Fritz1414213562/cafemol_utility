#include"../../include/misc/ContainerJudge.hpp"
#include<vector>
#include<array>
#include<iostream>

int main() {

	std::vector<std::array<int, 2>> v;
	v.push_back({0, 1});
	v.push_back({2, 5});
	v.push_back({3, 3});
	v.push_back({5, 1});
	v.push_back({9, 0});
	v.push_back({1, 2});

	std::vector<std::array<int, 2>> keys;
	keys.push_back({3, 3});
	keys.push_back({2, 2});
	keys.push_back({9, 0});
	keys.push_back({0, 9});
	keys.push_back({1, 2});
	keys.push_back({1, 2});
	for (const std::array<int, 2>& key : keys) {
		bool res = cafemol::library::is_containercontains(v, key);
//		std::cout << res << std::endl;
		if (res) {
			std::cout << key[0] << " " << key[1] << std::endl;
		}
	}

	return 0;
}

