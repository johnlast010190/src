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
    Foam::dynamicFreeSurfaceFvMesh

Description
    Mesh motion specifically for the free surface level of fluid tank.

Usage

dynamicFvMesh   dynamicFreeSurfaceFvMesh;
motionSolverLibs ( "libfvMotionSolvers.so" );
dynamicFreeSurfaceFvMeshCoeffs
{
    origin ( 0.5 1.5 0);     // origin
    direction (1 1 0);       // direction
    initialDisplacement 0.1; // initial compression height
    heightRef 0.3;           // undeformed region
    velocity constant 0.5;
}

SourceFiles
    dynamicFreeSurfaceFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicFreeSurfaceFvMesh_H
#define dynamicFreeSurfaceFvMesh_H

#include "dynamicFvMesh/dynamicFvMesh.H"
#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point/pointIOField.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicFreeSurfaceFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicFreeSurfaceFvMesh
:
    public dynamicFvMesh
{
    // Private data

        dictionary dynamicMeshCoeffs_;

    // Velocity
        autoPtr<Function1<scalar>> velocityProfile_;

    // Origin
       vector origin_;

    // Direction
       vector motionDir_;

    // Height below points are not moved
        scalar heightRef_;

    // Maximum mesh height along Dir
        scalar heightMax_;

    // Initial displacement
        scalar initialDisplacement_;

    //  Stationary points
        pointIOField stationaryPoints_;

    //  used to move points above reference height
        scalarField identifier_;

    //  used to ensure uniform mesh stretch
        scalarField ratio_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicFreeSurfaceFvMesh(const dynamicFreeSurfaceFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicFreeSurfaceFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicFreeSurfaceFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicFreeSurfaceFvMesh(const IOobject& io);


    //- Destructor
        ~dynamicFreeSurfaceFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
