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
    (c) 2021-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "strainRateLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "finiteVolume/fvc/fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(strainRateLaw, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        strainRateLaw,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::strainRateLaw::strainRateLaw
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    U_(nullptr),
    gradU_(nullptr)
{
    strainRateLaw::read();
}


Foam::autoPtr<Foam::strainRateLaw>
Foam::strainRateLaw::clone() const
{
    return autoPtr<strainRateLaw>
    (
        new strainRateLaw(*this)
    );
}


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

const volTensorField& Foam::strainRateLaw::gradU() const
{
    if (!gradU_.valid() || !gradU_().upToDate(*U_))
    {
        if (gradU_.valid() && gradU_().ownedByRegistry())
        {
            gradU_.release();
        }
        gradU_.reset(fvc::grad(*U_).ptr());
    }
    return *gradU_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::strainRateLaw::castScalarModel
(
    const word& modelName
)
{
    if (modelName == strainRateModel::typeName)
    {
        return dynamic_cast<strainRateModel*>(this);
    }
    return nullptr;
}


bool Foam::strainRateLaw::updateVectorField
(
    const word& fieldName,
    const volVectorField& volField
)
{
    if (fieldName == phasePropertyName("U", phaseName_))
    {
        U_ = &volField;
        return true;
    }
    return false;
}


Foam::scalar Foam::strainRateLaw::strainRateCell(const label celli) const
{
    return sqrt(2.0)*mag(symm(gradU()[celli]));
}


Foam::tmp<Foam::scalarField> Foam::strainRateLaw::strainRatePatch
(
    const label patchi
) const
{
    return sqrt(2.0)*mag(symm(gradU().boundaryField()[patchi]));
}


Foam::tmp<Foam::scalarField> Foam::strainRateLaw::strainRateInternal() const
{
    return sqrt(2.0)*mag(symm(gradU().primitiveField()));
}


Foam::tmp<Foam::volScalarField> Foam::strainRateLaw::strainRateGeometric
() const
{
    return sqrt(2.0)*mag(symm(gradU()));
}


bool Foam::strainRateLaw::read()
{
    U_ = lookupPtr<vector>("U");
    return true;
}


// ************************************************************************* //
