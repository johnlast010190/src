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
    Foam::geomDecomp

Description
    Geometrical domain decomposition

SourceFiles
    geomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef geomDecomp_H
#define geomDecomp_H

#include "decompositionMethod/decompositionMethod.H"
#include "primitives/Vector/Vector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class geomDecomp Declaration
\*---------------------------------------------------------------------------*/

class geomDecomp
:
    public decompositionMethod
{

protected:

    // Protected data

        const dictionary& geomDecomDict_;

        Vector<label> n_;
        scalar delta_;
        tensor rotDelta_;

public:

    // Constructors

        //- Construct given the decomposition dictionary
        //  and the derived type name
        geomDecomp
        (
            const dictionary& decompositionDict,
            const word& derivedType
        );

    // Member Functions

        Vector<label> n() const
        {
            return n_;
        }

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        ) = 0;

        //- Like decompose but with uniform weights on the points
        virtual labelList decompose(const pointField&) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
