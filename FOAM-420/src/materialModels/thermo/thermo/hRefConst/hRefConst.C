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

CopyrightCopyright
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hRefConst.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hRefConst, 0);
    addToRunTimeSelectionTable(materialModel, hRefConst, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hRefConst::hRefConst
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
    dep_.setSize(4);
    hRefConst::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hRefConst::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    const word tName = materialTables_.tableName(phaseName_, specieName_);

    // Model names
    const word& CpName = CpModel::typeName;
    const word& HaName = HaModel::typeName;
    const word& HsName = HsModel::typeName;
    const word& SName = SModel::typeName;
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
}


Foam::baseModels<Foam::scalar>* Foam::hRefConst::castScalarModel
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


Foam::scalar Foam::hRefConst::CpCell(const label celli) const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoCell(celli);
    }
    return
        CpThermoCell(celli) + sMod_[CpDeparture][celli];
}


Foam::tmp<Foam::scalarField> Foam::hRefConst::CpPatch
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


Foam::tmp<Foam::scalarField> Foam::hRefConst::CpInternal() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoInternal();
    }
    return CpThermoInternal() + sMod_[CpDeparture].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::hRefConst::CpGeometric() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoGeometric();
    }
    return CpThermoGeometric() + sMod_[CpDeparture]();
}


Foam::scalar Foam::hRefConst::HaCell(const label celli) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoCell(celli);
    }
    return
        sMod_[HDeparture][celli] + HaThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::hRefConst::HaPatch
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


Foam::tmp<Foam::scalarField> Foam::hRefConst::HaInternal() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoInternal();
    }
    return
        sMod_[HDeparture].primitiveField()
      + HaThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::hRefConst::HaGeometric() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoGeometric();
    }
    return
        sMod_[HDeparture]() + HaThermoGeometric();
}


Foam::scalar Foam::hRefConst::HsCell(const label celli) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HsThermoCell(celli);
    }
    return
        sMod_[HDeparture][celli] + HsThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::hRefConst::HsPatch
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


Foam::tmp<Foam::scalarField> Foam::hRefConst::HsInternal() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HsThermoInternal();
    }
    return
        sMod_[HDeparture].primitiveField()
      + HsThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::hRefConst::HsGeometric() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HsThermoGeometric();
    }
    return
        sMod_[HDeparture]() + HsThermoGeometric();
}


Foam::scalar Foam::hRefConst::SCell(const label celli) const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoCell(celli);
    }
    return
        sMod_[SDeparture][celli] + SThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::hRefConst::SPatch
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


Foam::tmp<Foam::scalarField> Foam::hRefConst::SInternal() const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoInternal();
    }
    return
        sMod_[SDeparture].primitiveField()
      + SThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::hRefConst::SGeometric() const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoGeometric();
    }
    return
        sMod_[SDeparture]() + SThermoGeometric();
}


bool Foam::hRefConst::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    Cp_ = dict_->lookup<scalar>("Cp");
    Hf_ = dict_->lookup<scalar>("Hf");
    Tref_ = dict_->lookup<scalar>("Tref");
    Href_ = dict_->lookup<scalar>("Href");

    return true;
}


// ************************************************************************* //
