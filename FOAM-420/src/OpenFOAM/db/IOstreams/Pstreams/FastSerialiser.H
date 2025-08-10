#ifndef FOAMCORE_FASTSERIALISER_H
#define FOAMCORE_FASTSERIALISER_H


#include "db/IOstreams/Sstreams/ISstream.H"
#include "db/IOstreams/Sstreams/OSstream.H"


namespace Foam {

template<typename Type>
class Field;

template<typename T>
class Map;

template<typename Key, typename Hash>
class HashSet;

class face;
class edge;


template<class Type1, class Type2>
class Tuple2;

/**
 * Serialise `thing` and write it to `os`.
 *
 * This serialisation mechanism is faster and more compact than any other serialisation mechanism FOAM seems
 * to currently have, but it is entirely non-portable. Don't expect to be able to save the
 * output and load it in a later version of the software, or on a machine with different endianness, float
 * representation, or word size.
 *
 * This is useful for transient serialisation needs, such as sending a message between processors during a
 * run. It is not appropriate to use this to serialise output data files, since there is (by design)
 * no support for format migrations.
 */
template<typename T>
void fastSerialise(Ostream& os, const T& thing);

template<typename T>
void fastDeserialise(Istream& is, T& thing);

template<typename T>
struct FastSerialiser {
    void serialise(Ostream& os, const T& thing) {
        if constexpr (contiguous<T>()) {
            // If it's memcpyable, we can just dump the bytes.
            os.write((const char*) &thing, sizeof(T));
        } else {
            // Otherwise, a specialisation must have been provided.
            static_assert(!std::is_same<T, T>::value, "Unimplemented fast serialiser. Write one, or make the affected type `contiguous`.");
        }
    }

    void deserialise(Istream& is, T& output) {
        if constexpr (contiguous<T>()) {
            is.read((char*) &output, sizeof(T));
        } else {
            static_assert(!std::is_same<T, T>::value, "Unimplemented fast deserialiser. Write one, or make the affected type `contiguous`.");
        }
    }
};

// The procedure is the same for any type of list, so let's have an alias party.
template<typename ListT>
void fastSerialiseList(Ostream& os, const ListT& thing) {
    label len = thing.size();
    fastSerialise(os, len);
    for (label i = 0; i < len; i++) {
        fastSerialise(os, thing[i]);
    }
}
template<typename ListT>
void fastDeserialiseList(Istream& is, ListT& output) {
    label len;
    fastDeserialise(is, len);
    output.setSize(len);
    for (label i = 0; i < len; i++) {
        fastDeserialise(is, output[i]);
    }
}

template<typename T>
struct FastSerialiser<List<T>> {
    void serialise(Ostream& os, const List<T>& thing) {
        fastSerialiseList(os, thing);
    }
    void deserialise(Istream& is, List<T>& output) {
        fastDeserialiseList(is, output);
    }
};

template<typename T, unsigned X, unsigned Y, unsigned Z>
struct FastSerialiser<DynamicList<T, X, Y, Z>> {
    void serialise(Ostream& os, const DynamicList<T, X, Y, Z>& thing) {
        fastSerialiseList(os, thing);
    }
    void deserialise(Istream& is, DynamicList<T, X, Y, Z>& output) {
        fastDeserialiseList(is, output);
    }
};

// FixedList doesn't need to store the length at all, of course.
template<typename T, unsigned int N>
struct FastSerialiser<FixedList<T, N>> {
    void serialise(Ostream& os, const FixedList<T, N>& thing) {
        for (label i = 0; i < (label) N; i++) {
            fastSerialise(os, thing[i]);
        }
    }
    void deserialise(Istream& is, FixedList<T, N>& output) {
        for (label i = 0; i < (label) N; i++) {
            fastDeserialise(is, output[i]);
        }
    }
};
template<typename T>
struct FastSerialiser<Field<T>> {
    void serialise(Ostream& os, const Field<T>& thing) {
        fastSerialiseList(os, thing);
    }
    void deserialise(Istream& is, Field<T>& output) {
        fastDeserialiseList(is, output);
    }
};
template<>
struct FastSerialiser<face> {
    void serialise(Ostream& os, const face& thing) {
        fastSerialiseList(os, (const List<label>&) thing);
    }
    void deserialise(Istream& is, face& output) {
        fastDeserialiseList(is, (List<label>&) output);
    }
};

template<>
struct FastSerialiser<edge> {
    void serialise(Ostream& os, const edge& thing) {
        fastSerialiseList(os, (const FixedList<label, 2>&) thing);
    }
    void deserialise(Istream& is, edge& output) {
        fastDeserialiseList(is, (FixedList<label, 2>&) output);
    }
};

// HashSets can be serialised just like lists. How much do we bet this will break
// things due to iteration order issues? :D
template<typename T, typename H>
struct FastSerialiser<HashSet<T, H>> {
    void serialise(Ostream& os, const HashSet<T, H>& thing) {
        label len = thing.size();
        fastSerialise(os, len);
        for (const auto it : thing) {
            fastSerialise(os, it.key());
        }
    }
    void deserialise(Istream& is, HashSet<T, H>& output) {
        label len;
        fastDeserialise(is, len);

        // TODO: Can probably be made more efficient using a batch-insertion routine on the map.
        for (label i = 0; i < len; i++) {
            T x; // TODO: Element uselessly copied. Need `emplace()`.
            fastDeserialise(is, x);
            output.insert(x);
        }
    }
};

template<typename V, typename K, typename H>
struct FastSerialiser<HashTable<V, K, H>> {
    void serialise(Ostream& os, const HashTable<V, K, H>& thing) {
        label len = thing.size();
        fastSerialise(os, len);
        for (auto it = thing.begin(); it != thing.end(); it++) {
            fastSerialise(os, it.key());
            fastSerialise(os, it.object());
        }
    }

    void deserialise(Istream& is, HashTable<V, K, H>& output) {
        label len;
        fastDeserialise(is, len);

        // TODO: Can probably be made more efficient using a batch-insertion routine on the map.
        for (label i = 0; i < len; i++) {
            K k; // TODO: Element uselessly copied. Need `emplace()`.
            V v; // TODO: Element uselessly copied. Need `emplace()`.
            fastDeserialise(is, k);
            fastDeserialise(is, v);
            output.insert(k, v);
        }
    }
};

// TODO: Ugh, make Map an alias already.
template<typename V>
struct FastSerialiser<Map<V>> {
    void serialise(Ostream& os, const Map<V>& thing) {
        label len = thing.size();
        fastSerialise(os, len);
        for (auto it = thing.begin(); it != thing.end(); it++) {
            fastSerialise(os, it.key());
            fastSerialise(os, it.object());
        }
    }

    void deserialise(Istream& is, Map<V>& output) {
        label len;
        fastDeserialise(is, len);

        // TODO: Can probably be made more efficient using a batch-insertion routine on the map.
        for (label i = 0; i < len; i++) {
            label k; // TODO: Element uselessly copied. Need `emplace()`.
            V v; // TODO: Element uselessly copied. Need `emplace()`.
            fastDeserialise(is, k);
            fastDeserialise(is, v);
            output.insert(k, v);
        }
    }
};

// Strings are easy, as are all the random string-like types that for some reason exist.
template<typename StringEsqueT>
void fastSerialiseString(Ostream& os, const StringEsqueT& thing) {
    label len = thing.size();
    fastSerialise(os, len);
    os.write(thing.data(), thing.size());
}

template<typename StringEsqueT>
void fastDeserialiseString(Istream& is, StringEsqueT& output) {
    label len;
    fastDeserialise(is, len);
    output.resize(len);
    is.read(output.data(), len);
}

// And now alises for all the things that are really just kinda stringy...
template<>
struct FastSerialiser<std::string> {
    void serialise(Ostream& os, const std::string& thing) {
        fastSerialiseString(os, thing);
    }
    void deserialise(Istream& is, std::string& output) {
        fastDeserialiseString(is, output);
    }
};
template<>
struct FastSerialiser<string> {
    void serialise(Ostream& os, const string& thing) {
        fastSerialiseString(os, thing);
    }
    void deserialise(Istream& is, string& output) {
        fastDeserialiseString(is, output);
    }
};
template<>
struct FastSerialiser<word> {
    void serialise(Ostream& os, const word& thing) {
        fastSerialiseString(os, thing);
    }
    void deserialise(Istream& is, word& output) {
        fastDeserialiseString(is, output);
    }
};
template<>
struct FastSerialiser<fileName> {
    void serialise(Ostream& os, const fileName& thing) {
        fastSerialiseString(os, thing);
    }
    void deserialise(Istream& is, fileName& output) {
        fastDeserialiseString(is, output);
    }
};

template<class Type1, class Type2>
struct FastSerialiser<Tuple2<Type1, Type2>> {
    void serialise(Ostream& os, const Tuple2<Type1, Type2>& thing) {
        fastSerialise(os, thing.first());
        fastSerialise(os, thing.second());
    }
    void deserialise(Istream& is, Tuple2<Type1, Type2>& output) {
        fastDeserialise(is, output.first());
        fastDeserialise(is, output.second());
    }
};

/**
 * Serialise `thing` and write it to `os`.
 *
 * This serialisation mechanism is faster and more compact than any other serialisation mechanism FOAM seems
 * to currently have, but it is entirely non-portable. Don't expect to be able to save the
 * output and load it in a later version of the software, or on a machine with different endianness, float
 * representation, or word size.
 *
 * This is useful for transient serialisation needs, such as sending a message between processors during a
 * run. It is not appropriate to use this to serialise output data files, since there is (by design)
 * no support for format migrations.
 */
template<typename T>
void fastSerialise(Ostream& os, const T& thing) {
    FastSerialiser<T>{}.serialise(os, thing);
}
template<typename T>
void fastDeserialise(Istream& is, T& thing) {
    FastSerialiser<T>{}.deserialise(is, thing);
}

} // namespace Foam


#endif //FOAMCORE_FASTSERIALISER_H
