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

#include "matSensibleEnthalpy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(matSensibleEnthalpy, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        matSensibleEnthalpy,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::matSensibleEnthalpy::matSensibleEnthalpy
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
    dep_.setSize(2);
}


Foam::autoPtr<Foam::matSensibleEnthalpy>
Foam::matSensibleEnthalpy::clone() const
{
    return autoPtr<matSensibleEnthalpy>(new matSensibleEnthalpy(*this));
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::matSensibleEnthalpy::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    if (modelName == CpvModel::typeName)
    {
        sMod_.set(Cp, models[CpModel::typeName]);
        dep_[0].model = models[CpvModel::typeName];
        dep_[0].dependencies.setSize(1);
        dep_[0].dependencies.set(0, models[CpModel::typeName]);
    }
    else if (modelName == HEModel::typeName)
    {
        sMod_.set(Hs, models[HsModel::typeName]);
        dep_[1].model = models[HEModel::typeName];
        dep_[1].dependencies.setSize(1);
        dep_[1].dependencies.set(0, models[HsModel::typeName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::matSensibleEnthalpy::castScalarModel
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


Foam::scalar Foam::matSensibleEnthalpy::CpvCell
(
    const label celli
) const
{
    return sMod_[Cp][celli];
}


Foam::tmp<Foam::scalarField> Foam::matSensibleEnthalpy::CpvPatch
(
    const label patchi
) const
{
    return sMod_[Cp].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::matSensibleEnthalpy::CpvInternal() const
{
    return sMod_[Cp].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::matSensibleEnthalpy::CpvGeometric() const
{
    return sMod_[Cp]();
}


Foam::scalar Foam::matSensibleEnthalpy::HECell
(
    const label celli
) const
{
    return sMod_[Hs][celli];
}


Foam::tmp<Foam::scalarField> Foam::matSensibleEnthalpy::HEPatch
(
    const label patchi
) const
{
    return sMod_[Hs].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::matSensibleEnthalpy::HEInternal() const
{
    return sMod_[Hs].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::matSensibleEnthalpy::HEGeometric() const
{
    return sMod_[Hs]();
}


// ************************************************************************* //
