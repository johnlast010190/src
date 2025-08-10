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
    Foam::fvMotionSolverEngineMesh

Description
    Foam::fvMotionSolverEngineMesh

SourceFiles
    fvMotionSolverEngineMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvMotionSolverEngineMesh_H
#define fvMotionSolverEngineMesh_H

#include "engineMesh/engineMesh/engineMesh.H"
#include "fvMotionSolvers/componentVelocity/componentLaplacian/velocityComponentLaplacianFvMotionSolver.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fvMotionSolverEngineMesh Declaration
\*---------------------------------------------------------------------------*/

class fvMotionSolverEngineMesh
:
    public engineMesh
{
    // Private data

        dimensionedScalar pistonLayers_;

        //- Mesh-motion solver to solve for the "z" component of the cell-centre
        //  displacements
        velocityComponentLaplacianFvMotionSolver motionSolver_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvMotionSolverEngineMesh(const fvMotionSolverEngineMesh&);

        //- Disallow default bitwise assignment
        void operator=(const fvMotionSolverEngineMesh&);


public:

    //- Runtime type information
    TypeName("fvMotionSolver");


    // Constructors

        //- Construct from IOobject
        fvMotionSolverEngineMesh(const IOobject& io);


    //- Destructor
    ~fvMotionSolverEngineMesh();


    // Member Functions

        // Edit

            virtual void moveMesh() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
