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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2020 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::dynamicMotionSolverFvMeshAMI

Description
    The dynamicMotionSolverFvMeshAMI

SourceFiles
    dynamicMotionSolverFvMeshAMI.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverFvMeshAMI_H
#define dynamicMotionSolverFvMeshAMI_H

#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                 Class dynamicMotionSolverFvMeshAMI Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverFvMeshAMI
:
    public dynamicFvMesh
{
    // Private data

        autoPtr<motionSolver> motionPtr_;


    // Private Member Functions

        //- No copy construct
        dynamicMotionSolverFvMeshAMI
        (
            const dynamicMotionSolverFvMeshAMI&
        ) = delete;

        //- No copy assignment
        void operator=(const dynamicMotionSolverFvMeshAMI&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverFvMeshAMI");


    // Constructors

        //- Construct from IOobject
        dynamicMotionSolverFvMeshAMI(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        dynamicMotionSolverFvMeshAMI
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );


    //- Destructor
    virtual ~dynamicMotionSolverFvMeshAMI() = default;


    // Member Functions

        //- Return the motionSolver
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
