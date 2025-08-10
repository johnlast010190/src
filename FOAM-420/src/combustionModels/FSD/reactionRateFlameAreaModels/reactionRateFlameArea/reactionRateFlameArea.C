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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "FSD/reactionRateFlameAreaModels/reactionRateFlameArea/reactionRateFlameArea.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reactionRateFlameArea, 0);
    defineRunTimeSelectionTable(reactionRateFlameArea, dictionary);
}

const Foam::fvMesh& Foam::reactionRateFlameArea::getMesh(const objectRegistry& obr) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return  dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return dynamic_cast<const fvMesh&>(obr);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRateFlameArea::reactionRateFlameArea
(
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr,
    const combustionModel& combModel
)
:
    coeffDict_(dict.optionalSubDict(modelType + "Coeffs")),
    obr_(obr),
    mesh_(getMesh(obr)),
    combModel_(combModel),
    fuel_(dict.lookup("fuel")),
    omega_
    (
        IOobject
        (
            "FSDomega",
            mesh_.time().timeName(),
            obr_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateFlameArea::~reactionRateFlameArea()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::reactionRateFlameArea::read(const dictionary& dict)
{
    dict.lookup("fuel") >> fuel_;

    return true;
}


// ************************************************************************* //
