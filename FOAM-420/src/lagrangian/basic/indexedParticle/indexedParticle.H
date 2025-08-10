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

Class
    Foam::indexedParticle

Description
    Adds label index to base particle

SourceFiles
    indexedParticle.H

\*---------------------------------------------------------------------------*/

#ifndef indexedParticle_H
#define indexedParticle_H

#include "particle/particle.H"
#include "db/IOstreams/IOstreams/IOstream.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class indexedParticle Declaration
\*---------------------------------------------------------------------------*/

class indexedParticle
:
    public particle
{
    // Private data

        label index_;


public:

    // Constructors

        //- Construct from components
        indexedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label index = 0
        )
        :
            particle(mesh, position, celli, tetFacei, tetPti),
            index_(index)
        {}

        //- Construct from components, with searching for tetFace and tetPt
        indexedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label index = 0
        )
        :
            particle(mesh, position, celli),
            index_(index)
        {}

        //- Construct from Istream
        indexedParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        )
        :
            particle(mesh, is, readFields)
        {}

        //- Construct as a copy
        indexedParticle(const indexedParticle& p)
        :
            particle(p)
        {}

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new indexedParticle(*this));
        }


    // Member functions

        label index() const
        {
            return index_;
        }

        label& index()
        {
            return index_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
