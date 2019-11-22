#ifndef ERROR_MESSAGE_HPP
#define ERROR_MESSAGE_HPP
#include<iostream>
#include<type_traits>
#include<string>
#include<tuple>

// dump to standard error output.

namespace cafemol {

namespace error_handling {

class Error_Output {

public:
	Error_Output() = default;
	~Error_Output() = default;

	template<typename ...Types>
	void output_Error(Types ...m_messages) {
		constexpr std::size_t Types_size = sizeof...(Types);
		static_assert(Types_size > 0, "No arguments");
		std::tuple<Types...> message_tuple = std::make_tuple(m_messages...);
		iterate_Message(message_tuple);

		std::exit(1);
	}

//	template<typename ...Args>
//	void output_Error(Args ...args) {
//		constexpr std::size_t Args_size = sizeof...(Args);
//		static_assert(Args_size > 0, "No arguments");
//		std::initializer_list<std::string> args_strings = {args...};	
//
//		bool is_initial_line = true;
//		for (const std::string& args_string : args_strings) {
//			if (is_initial_line) {
//				output_ErrorMessage(args_string);
//				is_initial_line = false;
//			}
//			else {
//				output_Message(args_string);
//			}
//		}
//
//		std::exit(1);
//	}
//
	template<typename ...Args>
	void operator()(Args ...args) {
		output_Error(args...);
	}

private:

	// add 'error_prefix' to message
	template<typename message_type>
	void output_ErrorMessage(const message_type& message) {
		std::cerr << error_prefix << message << std::endl;
	}
	// not add 'error_prefix'
	template<typename message_type>
	void output_Message(const message_type& message) {
		std::cerr << blank_prefix << message << std::endl;
	}

	template<std::size_t iterative_num = 0, typename Tuple_Type>
	void iterate_Message(const Tuple_Type& Tuple) {
		if constexpr ((iterative_num > 0) && (iterative_num < std::tuple_size<Tuple_Type>::value)) {
			output_Message(std::get<iterative_num>(Tuple));
			iterate_Message<iterative_num + 1>(Tuple);
		}
		else if constexpr (iterative_num == 0) {
			output_ErrorMessage(std::get<iterative_num>(Tuple));
			iterate_Message<iterative_num + 1>(Tuple);
		}
	}


//	// add 'error_prefix' to message
//	template<typename message_type>
//	void output_ErrorMessage(const message_type& message) {
//		std::cerr << error_prefix << message << std::endl;
//	}
//	// not add 'error_prefix'
//	template<typename message_type>
//	void output_Message(const message_type& message) {
//		std::cerr << blank_prefix << message << std::endl;
//	}
//
//
	const std::string error_prefix = "Error:  ";
	const std::string blank_prefix = "	";

};
}
}

#endif /* ERROR_MESSAGE_HPP */
