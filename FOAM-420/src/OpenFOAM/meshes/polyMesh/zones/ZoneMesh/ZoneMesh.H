/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016-2018 OpenCFD Ltd.

Class
    Foam::ZoneMesh

Description
    A list of mesh zones.

SourceFiles
    ZoneMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ZoneMesh_H
#define ZoneMesh_H

#include "containers/Lists/List/List.H"
#include "db/regIOobject/regIOobject.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "containers/HashTables/Map/Map.H"
#include "containers/Lists/PackedList/PackedBoolList.H"
#include "primitives/strings/wordRes/wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

template<class ZoneType, class MeshType> class ZoneMesh;

template<class ZoneType, class MeshType>
Ostream& operator<<(Ostream& os, const ZoneMesh<ZoneType, MeshType>& zones);

/*---------------------------------------------------------------------------*\
                           Class ZoneMesh Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType, class MeshType>
class ZoneMesh
:
    public PtrList<ZoneType>,
    public regIOobject
{
    // Private data

        //- Reference to mesh
        const MeshType& mesh_;

        //- Map of zone labels for given element
        mutable Map<label>* zoneMapPtr_;


    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool read();

        //- Create zone map
        void calcZoneMap() const;

        //- Templated implementation for names()
        template<class UnaryMatchPredicate>
        static wordList namesImpl
        (
            const PtrList<ZoneType>& zones,
            const UnaryMatchPredicate& matcher,
            const bool doSort
        );

        //- Templated implementation for indices()
        template<class UnaryMatchPredicate>
        static labelList indicesImpl
        (
            const PtrList<ZoneType>& list,
            const UnaryMatchPredicate& matcher
        );

        //- Disallow construct as copy
        ZoneMesh(const ZoneMesh&) = delete;

        //- Disallow assignment
        void operator=(const ZoneMesh<ZoneType, MeshType>&) = delete;


public:

    // Constructors

        //- Read constructor given IOobject and a MeshType reference
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh
        );

        //- Construct given size
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            const label size
        );

        //- Construct given a PtrList
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            const PtrList<ZoneType>& pzm
        );


    //- Destructor
    ~ZoneMesh();


    // Member Functions

        //- Return the mesh reference
        const MeshType& mesh() const
        {
            return mesh_;
        }

        //- Map of zones containing zone index for all zoned elements
        //  Return -1 if the object is not in the zone
        const Map<label>& zoneMap() const;

        //- Given a global object index, return the zone it is in.
        // If object does not belong to any zones, return -1
        label whichZone(const label objectIndex) const;

        //- Return a list of zone types
        wordList types() const;

        //- A list of the zone names
        wordList names() const;

        //- A list of the zone names satisfying the input matcher
        wordList names(const wordRe& matcher) const;

        //- A list of the zone names satisfying the input matchers
        wordList names(const wordRes& matcher) const;

        //- Sorted list of the zone names
        wordList sortedNames() const;

        //- Sorted list of the zone names satisfying the input matcher
        wordList sortedNames(const wordRe& matcher) const;

        //- Sorted list of the zone names satisfying the input matchers
        wordList sortedNames(const wordRes& matcher) const;


        //- Return zone indices for all matches
        labelList findIndices(const keyType& key) const;

        //- Return zone indices for all matches
        labelList findIndices(const wordRes& matcher) const;


        //- Return zone index for the first match, return -1 if not found
        label findIndex(const keyType& key) const;

        //- Find zone index given a name
        label findZoneID(const word& zoneName) const;

        //- Mark cells that match the zone specification
        PackedBoolList findMatching(const keyType& key) const;

        //- Return all elements (cells, faces, points) contained in the
        //- listed zones.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        PackedBoolList findMatching(const labelUList& zoneIds) const;

        //- Return all elements (cells, faces, points) that match the zone
        //- specification as a bitSet.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        PackedBoolList findMatching(const wordRes& matcher) const;


        //- Clear addressing
        void clearAddressing();

        //- Clear the zones
        void clear();

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have all zones and in same order. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const;

        //- Correct zone mesh after moving points
        void movePoints(const pointField& p);

        //- writeData member function required by regIOobject
        bool writeData(Ostream& os) const;

    // Member Operators

        //- Return const and non-const reference to ZoneType by index.
        using PtrList<ZoneType>::operator[];

        //- Return const reference to ZoneType by name.
        const ZoneType& operator[](const word& zoneName) const;

        //- Return reference to ZoneType by name.
        ZoneType& operator[](const word& zoneName);


    // Ostream operator

        friend Ostream& operator<< <ZoneType, MeshType>
        (
            Ostream& os,
            const ZoneMesh<ZoneType, MeshType>& zones
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshes/polyMesh/zones/ZoneMesh/ZoneMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
