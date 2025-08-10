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
    (c) 2010-2019 Esi Ltd.

Class
    Foam::mappedWallFvPatch

Description
    A region-coupled wall patch, which acts as an lduInterface for
    implicit coupling if used with a monolithicMatrix

SourceFiles
    mappedWallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedWallFvPatch_H
#define mappedWallFvPatch_H

#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "mappedPatches/mappedPolyPatch/mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mappedWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedWallFvPatch
:
    virtual public mappedFvPatch,
    virtual public wallFvPatch
{

protected:

    // Protected members

        virtual void initMovePoints()
        {
            wallFvPatch::initMovePoints();
            mappedFvPatch::initMovePoints();
        }

public:

    //- Runtime type information
    TypeName(mappedWallPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        mappedWallFvPatch
        (
            const polyPatch& pp,
            const fvBoundaryMesh& bm
        )
        :
            fvPatch(pp, bm),
            mappedFvPatch(pp, bm),
            wallFvPatch(pp, bm)
        {}


    //- Destructor
    virtual ~mappedWallFvPatch()
    {}


    // Member Functions

        // Access

            //- Return mappedWallFvPatch neighbour patch
            virtual const mappedWallFvPatch&
            nbrPatch() const
            {
                return refCast<const mappedWallFvPatch>
                (
                    nbrFvMesh().boundary()[nbrPatchID()]
                );
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
