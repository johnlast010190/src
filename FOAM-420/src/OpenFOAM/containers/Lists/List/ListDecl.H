#ifndef FOAMCORE_LISTDECL_H
#define FOAMCORE_LISTDECL_H

#include "primitives/ints/label/label.H"

// *Just* the declaration of List (and potentially some other forward declarations.
// List.H has become an insane tangle of nonsense including many definitions, meaning
// that many things that want to merely forward-declare List end up being included by it.
// To avoid repeating declarations over and over, we centralise them here.

namespace Foam
{


class Istream;
class Ostream;

template<class T> class List;
template<class T> class UList;
template<class T> Ostream& operator<<(Ostream&, const UList<T>&);
template<class T> Istream& operator>>(Istream&, UList<T>&);

class labelRange;
template<class T> class SubList;

using labelList = List<label>;
using labelUList = UList<label>;

// Forward declaration of friend functions and operators
template<class T> Istream& operator>>(Istream& is, List<T>& L);

template<class T, unsigned Size> class FixedList;
template<class T> class PtrList;

class SLListBase;
template<class LListBase, class T> class LList;
template<class T>
using SLList = LList<SLListBase, T>;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicList;

template<class T> class SortableList;
template<class T> class IndirectList;
template<class T> class UIndirectList;
template<class T> class BiIndirectList;

} // namespace Foam

#endif //FOAMCORE_LISTDECL_H
