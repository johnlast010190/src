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
    (c) 2016 OpenCFD ltd.

\*---------------------------------------------------------------------------*/

#include "maxDeltaxyzCubeRootLESDelta.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(maxDeltaxyzCubeRootLESDelta, 0);
    addToRunTimeSelectionTable
    (
        LESdelta,
        maxDeltaxyzCubeRootLESDelta,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyzCubeRootLESDelta::calcDelta()
{
    maxDeltaxyz_.calcDelta();
    cubeRootVolDelta_.calcDelta();

    delta_ =
        max
        (
            static_cast<const volScalarField&>(maxDeltaxyz_),
            static_cast<const volScalarField&>(cubeRootVolDelta_)
        );

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::maxDeltaxyzCubeRootLESDelta::maxDeltaxyzCubeRootLESDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    maxDeltaxyz_(name, turbulence, dict.subDict(typeName + "Coeffs")),
    cubeRootVolDelta_(name, turbulence, dict.subDict(typeName + "Coeffs"))
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyzCubeRootLESDelta::read(const dictionary& dict)
{
   maxDeltaxyz_.read(dict.subDict(typeName + "Coeffs"));
   cubeRootVolDelta_.read(dict.subDict(typeName + "Coeffs"));

   calcDelta();
}


void Foam::LESModels::maxDeltaxyzCubeRootLESDelta::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
