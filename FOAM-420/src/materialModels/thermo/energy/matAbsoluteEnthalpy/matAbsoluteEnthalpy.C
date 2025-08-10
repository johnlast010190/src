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

#include "matAbsoluteEnthalpy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(matAbsoluteEnthalpy, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        matAbsoluteEnthalpy,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::matAbsoluteEnthalpy::matAbsoluteEnthalpy
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


Foam::autoPtr<Foam::matAbsoluteEnthalpy>
Foam::matAbsoluteEnthalpy::clone() const
{
    return autoPtr<matAbsoluteEnthalpy>(new matAbsoluteEnthalpy(*this));
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::matAbsoluteEnthalpy::updateTable(const word& modelName)
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
        sMod_.set(Ha, models[HaModel::typeName]);
        dep_[1].model = models[HEModel::typeName];
        dep_[1].dependencies.setSize(1);
        dep_[1].dependencies.set(0, models[HaModel::typeName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::matAbsoluteEnthalpy::castScalarModel
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


Foam::scalar Foam::matAbsoluteEnthalpy::CpvCell
(
    const label celli
) const
{
    return sMod_[Cp][celli];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteEnthalpy::CpvPatch
(
    const label patchi
) const
{
    return sMod_[Cp].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteEnthalpy::CpvInternal() const
{
    return sMod_[Cp].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::matAbsoluteEnthalpy::CpvGeometric() const
{
    return sMod_[Cp]();
}


Foam::scalar Foam::matAbsoluteEnthalpy::HECell
(
    const label celli
) const
{
    return sMod_[Ha][celli];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteEnthalpy::HEPatch
(
    const label patchi
) const
{
    return sMod_[Ha].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::matAbsoluteEnthalpy::HEInternal() const
{
    return sMod_[Ha].primitiveField();
}


Foam::tmp<Foam::volScalarField>
Foam::matAbsoluteEnthalpy::HEGeometric() const
{
    return sMod_[Ha]();
}


// ************************************************************************* //
