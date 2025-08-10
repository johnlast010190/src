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
    (c) 2020 Esi Ltd.

Class
    Foam::multiSolidBodyGIBMotionSolver

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function, returning STL deformation for use in solidBodyGIBFvMesh
    motion solver.

SourceFiles
    multiSolidBodyGIBMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef multiSolidBodyGIBMotionSolver_H
#define multiSolidBodyGIBMotionSolver_H

#include "motionSolvers/motionSolver/motionSolver.H"
#include "solidBodyMotionFunctions/solidBodyMotionFunction/solidBodyMotionFunction.H"
#include "motionSolvers/displacement/solidBody/multiSolidBodyMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiSolidBodyGIBMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class multiSolidBodyGIBMotionSolver
:
    public motionSolver
{
    // Private data

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;

        //- name of the stl geometry
        const word triName_;

        //- Original stl points
        pointField triPoints0_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiSolidBodyGIBMotionSolver(const multiSolidBodyGIBMotionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const multiSolidBodyGIBMotionSolver&);

public:

    //- Runtime type information
    TypeName("multiSolidBodyGIBMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        multiSolidBodyGIBMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );


    //- Destructor
    ~multiSolidBodyGIBMotionSolver();


    // Member Functions

        //- Return point location of the STL obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&)
        {}

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&)
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
