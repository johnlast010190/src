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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "matAbsoluteInternalEnergy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(matAbsoluteInternalEnergy, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        matAbsoluteInternalEnergy,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::matAbsoluteInternalEnergy::matAbsoluteInternalEnergy
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(3);
}


Foam::autoPtr<Foam::matAbsoluteInternalEnergy>
Foam::matAbsoluteInternalEnergy::clone() const
{
    return autoPtr<matAbsoluteInternalEnergy>
    (
        new matAbsoluteInternalEnergy(*this)
    );
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::matAbsoluteInternalEnergy::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    if (modelName == CpvModel::typeName)
    {
        sMod_.set(Cv, models[CvModel::typeName]);
        dep_[0].model = models[CpvModel::typeName];
        dep_[0].dependencies.setSize(1);
        dep_[0].dependencies.set(0, models[CvModel::typeName]);
    }
    else if (modelName == CpByCpvModel::typeName)
    {
        sMod_.set(gamma, models[gammaModel::typeName]);
        dep_[1].model = models[CpByCpvModel::typeName];
        dep_[1].dependencies.setSize(1);
        dep_[1].dependencies.set(0, models[gammaModel::typeName]);
    }
    else if (modelName == HEModel::typeName)
    {
        sMod_.set(Ea, models[EaModel::typeName]);
        dep_[1].model = models[HEModel::typeName];
        dep_[1].dependencies.setSize(1);
        dep_[1].dependencies.set(0, models[EaModel::typeName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::matAbsoluteInternalEnergy::castScalarModel
(
    const word& modelName
)
{
    if (modelName == CpvModel::typeName)
    {
        return dynamic_cast<CpvModel*>(this);
    }
    else if (modelName == CpByCpvModel::typeName)
    {
        return dynamic_cast<CpByCpvModel*>(this);
    }
    else if (modelName == HEModel::typeName)
    {
        return dynamic_cast<HEModel*>(this);
    }
    return nullptr;
}


Foam::scalar Foam::matAbsoluteInternalEnergy::CpByCpvCell
(
    const label celli
) const
{
    return sMod_[gamma][celli];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteInternalEnergy::CpByCpvPatch
(
    const label patchi
) const
{
    return sMod_[gamma].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteInternalEnergy::CpByCpvInternal() const
{
    return sMod_[gamma].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::matAbsoluteInternalEnergy::CpByCpvGeometric() const
{
    return sMod_[gamma]();
}


Foam::scalar Foam::matAbsoluteInternalEnergy::CpvCell
(
    const label celli
) const
{
    return sMod_[Cv][celli];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteInternalEnergy::CpvPatch
(
    const label patchi
) const
{
    return sMod_[Cv].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteInternalEnergy::CpvInternal() const
{
    return sMod_[Cv].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::matAbsoluteInternalEnergy::CpvGeometric() const
{
    return sMod_[Cv]();
}


Foam::scalar Foam::matAbsoluteInternalEnergy::HECell
(
    const label celli
) const
{
    return sMod_[Ea][celli];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteInternalEnergy::HEPatch
(
    const label patchi
) const
{
    return sMod_[Ea].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteInternalEnergy::HEInternal() const
{
    return sMod_[Ea].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::matAbsoluteInternalEnergy::HEGeometric() const
{
    return sMod_[Ea]();
}


// ************************************************************************* //
