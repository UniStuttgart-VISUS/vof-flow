#pragma once

#include <type_traits>
#include <variant>

#include <boost/variant.hpp>

template<typename T, typename... Types>
inline std::add_pointer_t<const T> variant_get(const boost::variant<Types...>* var) {
    return boost::get<T>(var);
}

template<typename T, typename... Types>
inline std::add_pointer_t<const T> variant_get(const std::variant<Types...>* var) {
    return std::get_if<T>(var);
}

template<typename Visitor, typename... Types>
inline typename Visitor::result_type do_visit(Visitor&& vis, const boost::variant<Types...>& var) {
    return boost::apply_visitor(vis, var);
}

template<typename Visitor, typename... Types>
inline typename Visitor::result_type do_visit(Visitor&& vis, const std::variant<Types...>& var) {
    return std::visit(vis, var);
}
