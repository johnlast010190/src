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
    (c) 2013-2015 OpenFOAM Foundation

Class
    Foam::displacementMeshMoverMotionSolver

Description
    Mesh motion solver for an fvMesh. Uses externalDisplacementMeshMover
    to do the mesh motion.

SourceFiles
    displacementMeshMoverMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementMeshMoverMotionSolver_H
#define displacementMeshMoverMotionSolver_H

#include "motionSolvers/displacement/displacement/displacementMotionSolver.H"
#include "externalDisplacementMeshMover/externalDisplacementMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class displacementMeshMoverMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementMeshMoverMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        mutable autoPtr<externalDisplacementMeshMover> meshMoverPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        displacementMeshMoverMotionSolver
        (
            const displacementMeshMoverMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const displacementMeshMoverMotionSolver&);


public:

    //- Runtime type information
    TypeName("displacementMeshMover");


    // Constructors

        //- Construct from polyMesh and dictionary
        displacementMeshMoverMotionSolver(const polyMesh&, const dictionary&);


    //- Destructor
    ~displacementMeshMoverMotionSolver();


    // Member Functions

        externalDisplacementMeshMover& meshMover() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
