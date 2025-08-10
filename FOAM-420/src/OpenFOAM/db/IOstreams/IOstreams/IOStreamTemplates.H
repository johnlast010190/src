#ifndef FOAMCORE_IOSTREAMTEMPLATES_HPP
#define FOAMCORE_IOSTREAMTEMPLATES_HPP


#include <tuple>

// A place for IO operators that should be defined universally, and should be always inlineable
// for performance reasons. Generally pertains to STL entities.


namespace Foam
{


/*---------------------------------------------------------------------------*\
          Tuple IO
\*---------------------------------------------------------------------------*/


template<size_t... Is, typename... Ts>
void writeTuple(
    std::index_sequence<Is...>,
    Ostream& os,
    const std::tuple<Ts...>& tuple
) {
    ((os << std::get<Is>(tuple)), ...);
}

template<size_t... Is, typename... Ts>
void readTuple(
    std::index_sequence<Is...>,
    Istream& os,
    std::tuple<Ts...>& tuple
) {
    ((os >> std::get<Is>(tuple)), ...);
}


template<typename... Ts>
Ostream& operator<<(Ostream& os, const std::tuple<Ts...>& tuple) {
    writeTuple(std::index_sequence_for<Ts...>{}, os, tuple);
    return os;
}

template<typename... Ts>
Istream& operator>>(Istream& is, std::tuple<Ts...>& tuple) {
    readTuple(std::index_sequence_for<Ts...>{}, is, tuple);
    return is;
}

}


#endif //FOAMCORE_IOSTREAMTEMPLATES_HPP
