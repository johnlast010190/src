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
    Foam::cellQuality

Description
    Class calculates cell quality measures.

SourceFiles
    cellQuality.C

\*---------------------------------------------------------------------------*/

#ifndef cellQuality_H
#define cellQuality_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class cellQuality Declaration
\*---------------------------------------------------------------------------*/

class cellQuality
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellQuality(const cellQuality&);

        //- Disallow default bitwise assignment
        void operator=(const cellQuality&);


public:

    // Constructors

        //- Construct from mesh
        cellQuality(const polyMesh& mesh);


    // Destructor - default


    // Member Functions

        //- Return cell non-orthogonality
        tmp<scalarField> nonOrthogonality() const;

        //- Return cell skewness
        tmp<scalarField> skewness() const;

        //- Return face non-orthogonality
        tmp<scalarField> faceNonOrthogonality() const;

        //- Return face skewness
        tmp<scalarField> faceSkewness() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
