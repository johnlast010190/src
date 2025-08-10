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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::fvMotionSolver

Description
    Base class for fvMesh based motionSolvers.

SourceFiles
    fvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fvMotionSolver_H
#define fvMotionSolver_H

#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class velocityMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class fvMotionSolver
{

protected:

    // Protected data

        //- The fvMesh to be moved
        const fvMesh& fvMesh_;


    // Protected Member Functions

        //- Create the corresponding patch types for cellMotion from those
        //  of the given pointMotion
        template<class Type>
        wordList cellMotionBoundaryTypes
        (
            const typename GeometricField<Type, pointPatchField, pointMesh>::
            Boundary& pmUbf
        ) const;

public:

    //- Runtime type information
    ClassName("displacementMotionSolver");

    // Constructors

        //- Construct from polyMesh
        fvMotionSolver(const polyMesh&);


    // Member Functions

        //- Return reference to the fvMesh to be moved
        const fvMesh& mesh() const
        {
            return fvMesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMotionSolvers/fvMotionSolver/fvMotionSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
