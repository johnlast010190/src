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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::sixDoFRigidBodyGIBMotionSolver

Description
    sixDoF motion class for GIB. The class returns the new stl points based on
    the computed forces. Then the new points are used inside the GIB related
    fvMesh classes

SourceFiles
    sixDoFRigidBodyGIBMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyGIBMotionSolver_H
#define sixDoFRigidBodyGIBMotionSolver_H

#include "motionSolvers/motionSolver/motionSolver.H"
#include "sixDoFRigidBodyMotion/sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sixDoFRigidBodyGIBMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyGIBMotionSolver
:
    public motionSolver
{
    // Private data

        //- Six dof motion object
        autoPtr<sixDoFRigidBodyMotion> motionPtr_;

        wordReList patches_;

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

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- name of the stl geometry
        const word triName_;

        //- init stl points
        mutable pointField* initStlPointsPtr_;

        //- init stl points
        mutable pointField* stlPointsDisPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDoFRigidBodyGIBMotionSolver
        (
            const sixDoFRigidBodyGIBMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const sixDoFRigidBodyGIBMotionSolver&);

public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyGIBMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        sixDoFRigidBodyGIBMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~sixDoFRigidBodyGIBMotionSolver();


    // Member Functions

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&){};

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
