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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hPower.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hPower, 0);
    addToRunTimeSelectionTable(materialModel, hPower, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hPower::hPower
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
    dep_.setSize(5);
    hPower::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hPower::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    const word tName = materialTables_.tableName(phaseName_, specieName_);

    // Model names
    const word& CpName = CpModel::typeName;
    const word& HaName = HaModel::typeName;
    const word& HsName = HsModel::typeName;
    const word& SName = SModel::typeName;
    const word& HcName = HcModel::typeName;
    const word& HaThermoName = HaThermoModel::typeName;
    const word& CpDepartureName = CpDepartureModel::typeName;
    const word& HDepartureName = HDepartureModel::typeName;
    const word& SDepartureName = SDepartureModel::typeName;

    if (modelName == CpName && models.found(CpDepartureName))
    {
        sMod_.set(CpDeparture, models[CpDepartureName]);
        dep_[0].model = models[CpName];
        dep_[0].dependencies.setSize(1);
        dep_[0].dependencies.set(0, models[CpDepartureName]);
    }
    else if (modelName == HaName && models.found(HDepartureName))
    {
        sMod_.set(HDeparture, models[HDepartureName]);
        dep_[1].model = models[HaName];
        dep_[1].dependencies.setSize(1);
        dep_[1].dependencies.set(0, models[HDepartureName]);
    }
    else if (modelName == HsName && models.found(HDepartureName))
    {
        sMod_.set(HDeparture, models[HDepartureName]);
        dep_[2].model = models[HsName];
        dep_[2].dependencies.setSize(1);
        dep_[2].dependencies.set(0, models[HDepartureName]);
    }
    else if (modelName == SName && models.found(SDepartureName))
    {
        sMod_.set(SDeparture, models[SDepartureName]);
        dep_[3].model = models[SName];
        dep_[3].dependencies.setSize(1);
        dep_[3].dependencies.set(0, models[SDepartureName]);
    }
    else if (modelName == HaThermoName)
    {
        sMod_.set(HsInd, models[HsName]);
        sMod_.set(HcInd, models[HcName]);
        dep_[4].model = models[HaThermoName];
        dep_[4].dependencies.setSize(2);
        dep_[4].dependencies.set(0, models[HsName]);
        dep_[4].dependencies.set(0, models[HcName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::hPower::castScalarModel
(
    const word& modelName
)
{
    if (modelName == CpModel::typeName)
    {
        return dynamic_cast<CpModel*>(this);
    }
    else if (modelName == HaModel::typeName)
    {
        return dynamic_cast<HaModel*>(this);
    }
    else if (modelName == HsModel::typeName)
    {
        return dynamic_cast<HsModel*>(this);
    }
    else if (modelName == HcModel::typeName)
    {
        return dynamic_cast<HcModel*>(this);
    }
    else if (modelName == SModel::typeName)
    {
        return dynamic_cast<SModel*>(this);
    }
    return nullptr;
}


Foam::scalar Foam::hPower::CpCell(const label celli) const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoCell(celli);
    }
    return
        CpThermoCell(celli) + sMod_[CpDeparture][celli];
}


Foam::tmp<Foam::scalarField> Foam::hPower::CpPatch
(
    const label patchi
) const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoPatch(patchi);
    }
    return
        CpThermoPatch(patchi) + sMod_[CpDeparture].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::hPower::CpInternal() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoInternal();
    }
    return CpThermoInternal() + sMod_[CpDeparture].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::hPower::CpGeometric() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoGeometric();
    }
    return CpThermoGeometric() + sMod_[CpDeparture]();
}


Foam::scalar Foam::hPower::HaThermoCell(const label celli) const
{
    return sMod_[HsInd][celli] + sMod_[HcInd][celli];
}


Foam::tmp<Foam::scalarField> Foam::hPower::HaThermoPatch
(
    const label patchi
) const
{
    return
        sMod_[HsInd].boundaryField()[patchi]
      + sMod_[HcInd].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::hPower::HaThermoInternal() const
{
    return
        sMod_[HsInd].primitiveField()
      + sMod_[HcInd].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::hPower::HaThermoGeometric() const
{
    return sMod_[HsInd]() + sMod_[HcInd]();
}


Foam::scalar Foam::hPower::HaCell(const label celli) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoCell(celli);
    }
    return
        sMod_[HDeparture][celli] + HaThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::hPower::HaPatch
(
    const label patchi
) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoPatch(patchi);
    }
    return
        sMod_[HDeparture].boundaryField()[patchi]
      + HaThermoPatch(patchi);
}


Foam::tmp<Foam::scalarField> Foam::hPower::HaInternal() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoInternal();
    }
    return
        sMod_[HDeparture].primitiveField()
      + HaThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::hPower::HaGeometric() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoGeometric();
    }
    return
        sMod_[HDeparture]() + HaThermoGeometric();
}


Foam::scalar Foam::hPower::HsCell(const label celli) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HsThermoCell(celli);
    }
    return
        sMod_[HDeparture][celli] + HsThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::hPower::HsPatch
(
    const label patchi
) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HsThermoPatch(patchi);
    }
    return
        sMod_[HDeparture].boundaryField()[patchi]
      + HsThermoPatch(patchi);
}


Foam::tmp<Foam::scalarField> Foam::hPower::HsInternal() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HsThermoInternal();
    }
    return
        sMod_[HDeparture].primitiveField()
      + HsThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::hPower::HsGeometric() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HsThermoGeometric();
    }
    return
        sMod_[HDeparture]() + HsThermoGeometric();
}


Foam::scalar Foam::hPower::SCell(const label celli) const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoCell(celli);
    }
    return
        sMod_[SDeparture][celli] + SThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::hPower::SPatch
(
    const label patchi
) const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoPatch(patchi);
    }
    return
        sMod_[SDeparture].boundaryField()[patchi]
      + SThermoPatch(patchi);
}


Foam::tmp<Foam::scalarField> Foam::hPower::SInternal() const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoInternal();
    }
    return
        sMod_[SDeparture].primitiveField()
      + SThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::hPower::SGeometric() const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoGeometric();
    }
    return
        sMod_[SDeparture]() + SThermoGeometric();
}


bool Foam::hPower::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    c0_ = dict_->lookup<scalar>("C0");
    n0_ = dict_->lookup<scalar>("n0");
    Tref_ = dict_->lookup<scalar>("Tref");
    Hf_ = dict_->lookup<scalar>("Hf");

    return true;
}


// ************************************************************************* //
