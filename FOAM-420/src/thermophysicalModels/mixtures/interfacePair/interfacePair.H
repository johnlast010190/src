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
    (c) 2023 Esi Ltd.

Class
    Foam::interfacePair

SourceFiles
    interfacePair.C

\*---------------------------------------------------------------------------*/

#ifndef interfacePair_H
#define interfacePair_H

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class interfacePair Declaration
\*---------------------------------------------------------------------------*/

class interfacePair
:
    public Pair<word>
{
public:

    class hash
    :
        public Hash<interfacePair>
    {
    public:

        hash()
        {}

        label operator()(const interfacePair& key) const
        {
            return word::hash()(key.first()) + word::hash()(key.second());
        }
    };


    // Constructors

        //- Destructor
        interfacePair()
        {}

        //- Construct from the names directly
        interfacePair
        (
            const word& alpha1Name,
            const word& alpha2Name
        );

        //- Construct from the names taken from volumetric field
        interfacePair
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        );


    // Friend Operators

        //- Compare two interface pairs
        friend bool operator==
        (
            const interfacePair& a,
            const interfacePair& b
        )
        {
            return
            (
                ((a.first() == b.first()) && (a.second() == b.second()))
             || ((a.first() == b.second()) && (a.second() == b.first()))
            );
        }

        //- Compare two interface pairs
        friend bool operator!=
        (
            const interfacePair& a,
            const interfacePair& b
        )
        {
            return (!(a == b));
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
