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
    Foam::indirectConstraintFvPatch

Description
    Foam::indirectConstraintFvPatch

SourceFiles
    indirectConstraintFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef indirectConstraintFvPatch_H
#define indirectConstraintFvPatch_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "indirectConstraintPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class indirectConstraintFvPatch Declaration
\*---------------------------------------------------------------------------*/

class indirectConstraintFvPatch
:
    public fvPatch
{

public:

    //- Runtime type information
    TypeName(indirectConstraintPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        indirectConstraintFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm)
        {}

    // Member functions

        // Access

            //- Return true if this patch field is constrained
            virtual bool constrained() const
            {
                return true;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
