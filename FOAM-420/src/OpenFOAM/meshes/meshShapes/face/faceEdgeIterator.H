namespace Foam {


/**
 * Class template for walking the edges of anything "face-ish".
 *
 * Edges are traversed in point order. This is used to implement `walkEdges()` on face
 * classes, and should generally not be used directly (use a range-for loop on
 * `walkEdges()` instead).
 *
 * The various face implementations don't share a common base, although they share common
 * behaviour. Concepts are overrated, apparently ;)
 *
 * This gizmo can walk the edges of anything that provides `size()`, `faceEdge(label)`,
 * and is expected to successfully optimise nicely when the backing face's list-length
 * is a compile-time constant.
*/
template<typename FaceT>
class face_edge_walker {
    const FaceT& theFace;

    class edge_iterator {
        label i;
        const FaceT& theFace;

    public:
        explicit edge_iterator(const FaceT& theFace, label i = 0):
                i(i), theFace(theFace) {}

        edge operator*() const {
            return theFace.faceEdge(i);
        }

        bool operator==(const edge_iterator& other) const;
        bool operator!=(const edge_iterator& other) const;

        edge_iterator& operator++() {
            i++;
            return *this;
        }
    };

public:
    explicit face_edge_walker(const FaceT& theFace):
            theFace(theFace) {}

    edge_iterator begin() const {
        return edge_iterator{theFace, 0};
    }
    edge_iterator end() const {
        return edge_iterator{theFace, theFace.size()};
    }
};


} // namespace Foam

template<typename T>
bool Foam::face_edge_walker<T>::edge_iterator::operator==(const Foam::face_edge_walker<T>::edge_iterator& other) const {
    return i == other.i;
}

template<typename T>
bool Foam::face_edge_walker<T>::edge_iterator::operator!=(const Foam::face_edge_walker<T>::edge_iterator& other) const {
    return i != other.i;
}
