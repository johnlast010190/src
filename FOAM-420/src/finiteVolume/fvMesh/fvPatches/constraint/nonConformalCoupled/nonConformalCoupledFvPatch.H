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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::nonConformalCoupledFvPatch

Description
    Non-conformal coupled FV patch. As nonConformalFvPatch, but is also
    coupled to another non-conformal patch.

See also
    Foam::nonConformalCoupledPolyPatch
    Foam::nonConformalFvPatch

SourceFiles
    nonConformalCoupledFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCoupledFvPatch_H
#define nonConformalCoupledFvPatch_H

#include "fvMesh/fvPatches/basic/coupled/coupledFvPatch.H"
#include "nonConformal/polyPatches/nonConformalCoupled/nonConformalCoupledPolyPatch.H"
#include "fvMesh/fvPatches/constraint/nonConformal/nonConformalFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalErrorFvPatch;

/*---------------------------------------------------------------------------*\
                  Class nonConformalCoupledFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCoupledFvPatch
:
    public nonConformalFvPatch
{
    // Private data

        //- Reference to the fvPatch
        const coupledFvPatch& patch_;

        //- Reference to the polyPatch
        const nonConformalCoupledPolyPatch& nonConformalCoupledPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights
        (
            scalarField& w,
            const vectorField& nbrSf,
            const vectorField& nbrDelta
        ) const;


public:

    //- Runtime type information
    TypeName(nonConformalCoupledPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch
        nonConformalCoupledFvPatch(const fvPatch& patch);


    //- Destructor
    virtual ~nonConformalCoupledFvPatch();


    // Member Functions

        // Access

            //- Return the reference to the polyPatch
            const nonConformalCoupledPolyPatch&
                nonConformalCoupledPatch() const
            {
                return nonConformalCoupledPolyPatch_;
            }

            //- Does this side own the patch?
            bool owner() const
            {
                return nonConformalCoupledPolyPatch_.owner();
            }

            //- Does the coupled side own the patch?
            bool neighbour() const
            {
                return nonConformalCoupledPolyPatch_.neighbour();
            }

            //- Return transformation between the coupled patches
            const transformer& transform() const
            {
                return nonConformalCoupledPolyPatch_.transform();
            }

            //- Error patch name
            const word& errorPatchName() const
            {
                return nonConformalCoupledPolyPatch_.errorPatchName();
            }

            //- Error patch ID
            label errorPatchID() const
            {
                return nonConformalCoupledPolyPatch_.errorPatchID();
            }

            //- Error patch
            const nonConformalErrorFvPatch& errorPatch() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
