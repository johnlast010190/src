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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::geometricSurfacePatch

Description
    The geometricSurfacePatch is like patchIdentifier but for surfaces.
    Holds type, name and index.

SourceFiles
    geometricSurfacePatch.C

\*---------------------------------------------------------------------------*/

#ifndef geometricSurfacePatch_H
#define geometricSurfacePatch_H

#include "primitives/strings/word/word.H"
#include "primitives/ints/label/label.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators
class geometricSurfacePatch;

Istream& operator>>(Istream& is, geometricSurfacePatch& p);
Ostream& operator<<(Ostream& os, const geometricSurfacePatch& p);

/*---------------------------------------------------------------------------*\
                    Class geometricSurfacePatch Declaration
\*---------------------------------------------------------------------------*/

class geometricSurfacePatch
{
    // Private data

        //- Type name of patch
        word geometricType_;

        //- Name of patch
        word name_;

        //- Index of patch in boundary
        label index_;

public:

    // Public data

        //- The name for an 'empty' type
        static const word emptyType;


    //- Runtime type information
    ClassName("geometricSurfacePatch");


    // Constructors

        //- Construct null
        geometricSurfacePatch();

        //- Construct null with specified index
        explicit geometricSurfacePatch(const label index);

        //- Construct from components
        geometricSurfacePatch
        (
            const word& name,
            const label index,
            const word& geometricType = word::null
        );

        //- Construct from components
        geometricSurfacePatch
        (
            const word& geometricType,
            const word& name,
            const label index
        );

        //- Construct from dictionary
        geometricSurfacePatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Construct from Istream
        geometricSurfacePatch(Istream& is, const label index);


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return name
        word& name()
        {
            return name_;
        }

        //- Return the geometric type of the patch
        const word& geometricType() const
        {
            return geometricType_;
        }

        //- Return the geometric type of the patch for modification
        word& geometricType()
        {
            return geometricType_;
        }

        //- Return the index of this patch in the surface mesh
        label index() const
        {
            return index_;
        }

        //- Return the index of this patch in the surface mesh for modification
        label& index()
        {
            return index_;
        }

        //- Write
        void write(Ostream& os) const;

        //- Write dictionary
        void writeDict(Ostream& os) const;


    // Ostream Operator

        friend Istream& operator>>(Istream& is, geometricSurfacePatch& p);
        friend Ostream& operator<<(Ostream& os, const geometricSurfacePatch& p);
};


// Global Operators

//- Compare patches for equality
bool operator==(const geometricSurfacePatch& a, const geometricSurfacePatch& b);

//- Compare patches for inequality
bool operator!=(const geometricSurfacePatch& a, const geometricSurfacePatch& b);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// Sadface.
#include "db/IOstreams/Pstreams/FastSerialiser.H"

namespace Foam {

template<>
struct FastSerialiser<geometricSurfacePatch> {
    void serialise(Ostream& os, const geometricSurfacePatch& thing) {
        fastSerialise(os, thing.name());
        fastSerialise(os, thing.geometricType());
    }
    void deserialise(Istream& is, geometricSurfacePatch& output) {
        fastDeserialise(is, output.name());
        fastDeserialise(is, output.geometricType());
    }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
