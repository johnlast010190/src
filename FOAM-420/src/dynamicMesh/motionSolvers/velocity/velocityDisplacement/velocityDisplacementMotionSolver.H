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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::velocityDisplacementMotionSolver

Description
    Mesh motion solver for a polyMesh. Wraps a displacement motion solver in a
    velocity motion solver.

SourceFiles
    velocityDisplacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityDisplacementMotionSolver_H
#define velocityDisplacementMotionSolver_H

#include "motionSolvers/velocity/velocityMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class displacementMotionSolver;

/*---------------------------------------------------------------------------*\
              Class velocityDisplacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class velocityDisplacementMotionSolver
:
    public velocityMotionSolver
{
    // Private data

        //- Displacement motion solver
        autoPtr<displacementMotionSolver> displacementMotionSolverPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        velocityDisplacementMotionSolver
        (
            const velocityDisplacementMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const velocityDisplacementMotionSolver&);

        //- Get the boundary condition types for the point displacement
        wordList pointDisplacementBoundaryTypes() const;


public:

    //- Runtime type information
    TypeName("velocityDisplacement");


    // Constructors

        //- Construct from polyMesh and dictionary
        velocityDisplacementMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );


    //- Destructor
    ~velocityDisplacementMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update geometry
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
