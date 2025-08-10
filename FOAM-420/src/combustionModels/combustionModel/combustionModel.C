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

#include "combustionModel/combustionModel.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combustionModel, 0);
}

const Foam::word Foam::combustionModel::combustionPropertiesName
(
    "combustionProperties"
);

const Foam::fvMesh& Foam::combustionModel::getMesh(const objectRegistry& obr) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return  dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return dynamic_cast<const fvMesh&>(obr);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModel::combustionModel
(
    const word& modelType,
    const objectRegistry& obr,
    const word& combustionProperties,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(combustionProperties, phaseName),
            getMesh(obr).time().constant(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    turbulencePtr_(),
    obr_(obr),
    mesh_(getMesh(obr)),
    active_(lookupOrDefault<Switch>("active", true)),
    coeffs_(optionalSubDict(modelType + "Coeffs")),
    modelType_(modelType),
    phaseName_(phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModel::~combustionModel()
{
    if (turbulencePtr_)
    {
        turbulencePtr_ = 0;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::combustionModel::read()
{
    if (regIOobject::read())
    {
        this->lookup("active") >> active_;
        coeffs_ = optionalSubDict(modelType_ + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
