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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "janaf.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(janaf, 0);
    addToRunTimeSelectionTable(materialModel, janaf, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::janaf::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorInFunction
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::janaf::janaf
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    isConverted_(false)
{
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(4);
    janaf::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::janaf::updateTable(const word& modelName)
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
    const word& CpDepartureName = CpDepartureModel::typeName;
    const word& HDepartureName = HDepartureModel::typeName;
    const word& SDepartureName = SDepartureModel::typeName;
    const word& RName = RModel::typeName;

    if (modelName == CpName && models.found(CpDepartureName))
    {
        sMod_.set(CpDeparture, models[CpDepartureName]);
        dep_[0].model = models[CpName];
        dep_[0].dependencies.setSize(2);
        dep_[0].dependencies.set(0, models[CpDepartureName]);
        dep_[0].dependencies.set(1, models[RName]);
    }
    else if (modelName == HaName && models.found(HDepartureName))
    {
        sMod_.set(HDeparture, models[HDepartureName]);
        dep_[1].model = models[HaName];
        dep_[1].dependencies.setSize(2);
        dep_[1].dependencies.set(0, models[HDepartureName]);
        dep_[1].dependencies.set(1, models[RName]);
    }
    else if (modelName == HsName)
    {
        sMod_.set(HaAll, models[HaName]);
        sMod_.set(HcAll, models[HcName]);
        dep_[2].model = models[HsName];
        dep_[2].dependencies.setSize(3);
        dep_[2].dependencies.set(0, models[HaName]);
        dep_[2].dependencies.set(1, models[HcName]);
        dep_[2].dependencies.set(2, models[RName]);
    }
    else if (modelName == SName && models.found(SDepartureName))
    {
        sMod_.set(SDeparture, models[SDepartureName]);
        dep_[3].model = models[SName];
        dep_[3].dependencies.setSize(1);
        dep_[3].dependencies.set(0, models[SDepartureName]);
        dep_[3].dependencies.set(1, models[RName]);
    }

    //- Needs to be set always
    sMod_.set(R, models[RName]);
    R_ = sMod_[R][0];

    if (!isConverted_)
    {
        isConverted_ = true;
        // Convert coefficients to mass-basis
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] *= R_;
            lowCpCoeffs_[coefLabel] *= R_;
        }
        checkInputData();
    }
}


Foam::baseModels<Foam::scalar>* Foam::janaf::castScalarModel
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


Foam::scalar Foam::janaf::CpCell(const label celli) const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoCell(celli);
    }
    return
        CpThermoCell(celli) + sMod_[CpDeparture][celli];
}


Foam::tmp<Foam::scalarField> Foam::janaf::CpPatch
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


Foam::tmp<Foam::scalarField> Foam::janaf::CpInternal() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoInternal();
    }
    return CpThermoInternal() + sMod_[CpDeparture].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::janaf::CpGeometric() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return CpThermoGeometric();
    }
    return CpThermoGeometric() + sMod_[CpDeparture]();
}


Foam::scalar Foam::janaf::HaCell(const label celli) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoCell(celli);
    }
    return
        sMod_[HDeparture][celli] + HaThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::janaf::HaPatch
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


Foam::tmp<Foam::scalarField> Foam::janaf::HaInternal() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoInternal();
    }
    return
        sMod_[HDeparture].primitiveField()
      + HaThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::janaf::HaGeometric() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return HaThermoGeometric();
    }
    return
        sMod_[HDeparture]() + HaThermoGeometric();
}


Foam::scalar Foam::janaf::HsCell(const label celli) const
{
    return sMod_[HaAll][celli] - sMod_[HcAll][celli];
}


Foam::tmp<Foam::scalarField> Foam::janaf::HsPatch
(
    const label patchi
) const
{
    return
        sMod_[HaAll].boundaryField()[patchi]
      - sMod_[HcAll].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::janaf::HsInternal() const
{
    return
        sMod_[HaAll].primitiveField()
      - sMod_[HcAll].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::janaf::HsGeometric() const
{
    return sMod_[HaAll]() - sMod_[HcAll]();
}


Foam::scalar Foam::janaf::SCell(const label celli) const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoCell(celli);
    }
    return
        sMod_[SDeparture][celli] + SThermoCell(celli);
}


Foam::tmp<Foam::scalarField> Foam::janaf::SPatch
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


Foam::tmp<Foam::scalarField> Foam::janaf::SInternal() const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoInternal();
    }
    return
        sMod_[SDeparture].primitiveField()
      + SThermoInternal();
}


Foam::tmp<Foam::volScalarField> Foam::janaf::SGeometric() const
{
    if (sMod_(SDeparture) == nullptr)
    {
        return SThermoGeometric();
    }
    return
        sMod_[SDeparture]() + SThermoGeometric();
}


bool Foam::janaf::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    Tlow_ = dict_->lookup<scalar>("Tlow");
    Thigh_ = dict_->lookup<scalar>("Thigh");
    Tcommon_ = dict_->lookup<scalar>("Tcommon");
    highCpCoeffs_ = dict_->lookup<scalarField>("highCpCoeffs");
    lowCpCoeffs_ = dict_->lookup<scalarField>("lowCpCoeffs");
    if (sMod_(R) != nullptr)
    {
        sMod_[R].read();
        R_ = sMod_[R][0];
    }

    return true;
}


// ************************************************************************* //
