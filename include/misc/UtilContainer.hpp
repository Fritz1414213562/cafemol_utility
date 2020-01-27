#ifndef UTIL_CONTAINER_HPP
#define UTIL_CONTAINER_HPP

namespace cafemol {

template<typename Type>
struct ValueWithIndex {

public:
	ValueWithIndex(const Type& input_index, const Type& input_value) :
		index(input_index), value(input_value) {}

	const Type& Index() {return index;}
	const Type& Value() {return value;}

//	Type Index() const {return index;}
//	Type Value() const {return value;}

private:
	const Type index;
	const Type value;

//	Type index;
//	Type value;
};


template<typename Type>
bool compare_VWI_Index(const ValueWithIndex<Type>& lhs, const ValueWithIndex<Type>& rhs) {
	return lhs.Index() < rhs.Index();
}

template<typename Type>
bool compare_VWI_Value(const ValueWithIndex<Type>& lhs, const ValueWithIndex<Type>& rhs) {
	return lhs.Value() < rhs.Value();
}


}

#endif /* UTIL_CONTAINER_HPP */
