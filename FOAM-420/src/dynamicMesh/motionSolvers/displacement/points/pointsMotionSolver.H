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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2023 Esi Ltd.

Class
    Foam::pointsMotionSolver

Description
    Virtual base class for displacement motion solvers, where the point
    motion is relative to a set of fixed points (points0).

SourceFiles
    pointsMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pointsMotionSolver_H
#define pointsMotionSolver_H

#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "meshes/primitiveShapes/point/pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                     Class pointsMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class pointsMotionSolver
:
    public motionSolver
{
protected:

    // Protected data

        //- Starting points
        pointIOField points0_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointsMotionSolver(const pointsMotionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const pointsMotionSolver&);


public:

    //- Runtime type information
    TypeName("pointsMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        pointsMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const word& type
        );

        //- Construct from mesh, dictionary and fixed points field
        pointsMotionSolver
        (
            const polyMesh&,
            const dictionary&,
            const pointIOField& points0,
            const word& type
        );


    //- Destructor
    virtual ~pointsMotionSolver();


    // Member Functions

        //- Return IO object for points0
        static IOobject points0IO(const polyMesh& mesh);

        //- Return reference to the reference field
        pointField& points0()
        {
            return points0_;
        }

        //- Return reference to the reference field
        const pointField& points0() const
        {
            return points0_;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&)
        {}

        //- Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&)
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
