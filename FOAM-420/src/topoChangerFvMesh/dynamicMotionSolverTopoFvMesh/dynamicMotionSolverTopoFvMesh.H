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
    (c) 2016 OpenCFD Ltd

Class
    Foam::dynamicMotionSolverTopoFvMesh

Description
    Sample topoChangerFvMesh that moves an object using a motion solver

SourceFiles
    dynamicMotionSolverTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverTopoFvMesh_H
#define dynamicMotionSolverTopoFvMesh_H

#include "topoChangerFvMesh/topoChangerFvMesh.H"
#include "motionSolvers/motionSolver/motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                      Class dynamicMotionSolverTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverTopoFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;

        //- Motion solver
        autoPtr<motionSolver> motionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicMotionSolverTopoFvMesh(const dynamicMotionSolverTopoFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotionSolverTopoFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverTopoFvMesh");


    // Constructors

        //- Construct from database
        explicit dynamicMotionSolverTopoFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicMotionSolverTopoFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
