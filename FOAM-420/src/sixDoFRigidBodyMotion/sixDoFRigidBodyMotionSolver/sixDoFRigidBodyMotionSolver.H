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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::sixDoFRigidBodyMotionSolver

Group
   grpMeshMotionSolvers

Description
    6-DoF solid-body mesh motion solver for an fvMesh.

    Applies SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    sixDoFRigidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionSolver_H
#define sixDoFRigidBodyMotionSolver_H

#include "motionSolvers/displacement/displacement/displacementMotionSolver.H"
#include "sixDoFRigidBodyMotion/sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sixDoFRigidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Six dof motion object
        sixDoFRigidBodyMotion motion_;

        wordReList patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Switch for test-mode in which only the
        //  gravitational body-force is applied
        Switch test_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDoFRigidBodyMotionSolver
        (
            const sixDoFRigidBodyMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const sixDoFRigidBodyMotionSolver&);


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        sixDoFRigidBodyMotionSolver
        (
            const polyMesh&,
            const dictionary& dict
        );


    //- Destructor
    ~sixDoFRigidBodyMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Write motion state information for restart
        virtual bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
