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
    (c) 2022 Esi Ltd.

Class
    Foam::cyclicAMIPointPatch

Description
    Cyclic AMI point patch - place holder only

SourceFiles
    cyclicAMIPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIPointPatch_H
#define cyclicAMIPointPatch_H

#include "meshes/pointMesh/pointPatches/derived/coupled/coupledFacePointPatch.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"
#include "meshes/pointMesh/pointBoundaryMesh/pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicAMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic AMI patch
        const cyclicAMIPolyPatch& cyclicAMIPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        cyclicAMIPointPatch(const cyclicAMIPointPatch&);

        //- Disallow default assignment
        void operator=(const cyclicAMIPointPatch&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName(cyclicAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicAMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~cyclicAMIPointPatch();


    // Member Functions

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return type();
        }

        //- Return the underlying cyclicAMIPolyPatch
        const cyclicAMIPolyPatch& cyclicAMIPatch() const
        {
            return cyclicAMIPolyPatch_;
        }

        //- Return neighbour point patch
        const cyclicAMIPointPatch& nbrPatch() const
        {
            label patchi = cyclicAMIPolyPatch_.nbrPatchID();
            const pointPatch& pp = this->boundaryMesh()[patchi];
            return refCast<const cyclicAMIPointPatch>(pp);
        }

        //- Return transformation between the coupled patches
        virtual const transformer& transform() const
        {
            return cyclicAMIPolyPatch_.transform();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
