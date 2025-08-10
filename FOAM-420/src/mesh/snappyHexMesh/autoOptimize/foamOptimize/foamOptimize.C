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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "autoOptimize/foamOptimize/foamOptimize.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(foamOptimize, 0);
    addToRunTimeSelectionTable
    (
        autoOptimize,
        foamOptimize,
        dictionary
    );
}


void Foam::foamOptimize::movePoints(Foam::pointField& newPoints)
{
    optimPtr_.reset
    (
        new hessianMeshOptimization
        (mesh_,newPoints, coeffsDict_, true)
    );
    newPoints = optimPtr_->newPoints();
}


void Foam::foamOptimize::optimize()
{
    optimPtr_.reset
    (
        new hessianMeshOptimization
        (mesh_, coeffsDict_, true)
    );
    tmp<pointField> newPoints = optimPtr_->newPoints();

    mesh_.movePoints(newPoints);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::foamOptimize::foamOptimize
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    autoOptimize(mesh, dict),
    mesh_(mesh),
    coeffsDict_(dict.subDict("foamOptimizeCoeffs")),
    optimPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamOptimize::~foamOptimize()
{}

// ************************************************************************* //
