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

#include "eConst.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eConst, 0);
    addToRunTimeSelectionTable(materialModel, eConst, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eConst::eConst
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name),
    Cv_("Cv", dimEnergy/(dimMass*dimTemperature), 0.0),
    Hf_("Hf", dimEnergy/dimMass, 0.0)
{
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(4);
    eConst::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eConst::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    const word tName = materialTables_.tableName(phaseName_, specieName_);

    // Model names
    const word& CpName = CpModel::typeName;
    const word& HaName = HaModel::typeName;
    const word& HsName = HsModel::typeName;
    const word& SName = SModel::typeName;
    const word& CpMCvName = CpMCvModel::typeName;
    const word& CpDepartureName = CpDepartureModel::typeName;
    const word& HDepartureName = HDepartureModel::typeName;
    const word& SDepartureName = SDepartureModel::typeName;

    if (modelName == CpName)
    {
        sMod_.set(CpMCv, models[CpMCvName]);
        dep_[0].model = models[CpName];
        dep_[0].dependencies.setSize(1);
        dep_[0].dependencies.set(0, models[CpMCvName]);

        if (models.found(CpDepartureName))
        {
            sMod_.set(CpDeparture, models[CpDepartureName]);
            dep_[0].dependencies.setSize(2);
            dep_[0].dependencies.set(0, models[CpMCvName]);
            dep_[0].dependencies.set(1, models[CpDepartureName]);
        }
    }
    else if (modelName == HaName)
    {
        sMod_.set(Cp, models[CpName]);
        dep_[1].model = models[HaName];
        dep_[1].dependencies.setSize(1);
        dep_[1].dependencies.set(0, models[CpName]);

        if (models.found(HDepartureName))
        {
            sMod_.set(HDeparture, models[HDepartureName]);
            dep_[1].dependencies.setSize(2);
            dep_[1].dependencies.set(0, models[CpName]);
            dep_[1].dependencies.set(1, models[HDepartureName]);
        }
    }
    else if (modelName == HsName)
    {
        sMod_.set(Cp, models[CpName]);
        dep_[2].model = models[HsName];
        dep_[2].dependencies.setSize(1);
        dep_[2].dependencies.set(0, models[CpName]);

        if (models.found(HDepartureName))
        {
            sMod_.set(HDeparture, models[HDepartureName]);
            dep_[2].dependencies.setSize(2);
            dep_[2].dependencies.set(0, models[CpName]);
            dep_[2].dependencies.set(1, models[HDepartureName]);
        }
    }
    else if (modelName == SName)
    {
        sMod_.set(Cp, models[CpName]);
        dep_[3].model = models[SName];
        dep_[3].dependencies.setSize(1);
        dep_[3].dependencies.set(0, models[CpName]);

        if (models.found(SDepartureName))
        {
            sMod_.set(SDeparture, models[SDepartureName]);
            dep_[3].dependencies.setSize(2);
            dep_[3].dependencies.set(0, models[CpName]);
            dep_[3].dependencies.set(1, models[SDepartureName]);
        }
    }
}


Foam::baseModels<Foam::scalar>* Foam::eConst::castScalarModel
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


Foam::scalar Foam::eConst::CpCell(const label celli) const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return Cv_.value() + sMod_[CpMCv][celli];
    }
    return
        Cv_.value() + sMod_[CpMCv][celli] + sMod_[CpDeparture][celli];
}


Foam::tmp<Foam::scalarField> Foam::eConst::CpPatch
(
    const label patchi
) const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return Cv_.value() + sMod_[CpMCv].boundaryField()[patchi];
    }
    return
         Cv_.value()
       + sMod_[CpMCv].boundaryField()[patchi]
       + sMod_[CpDeparture].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::eConst::CpInternal() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return Cv_.value() + sMod_[CpMCv].primitiveField();
    }
    return
        Cv_.value()
      + sMod_[CpMCv].primitiveField()
      + sMod_[CpDeparture].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::eConst::CpGeometric() const
{
    if (sMod_(CpDeparture) == nullptr)
    {
        return Cv_ + sMod_[CpMCv]();
    }
    return Cv_ + sMod_[CpMCv]() + sMod_[CpDeparture]();
}


Foam::scalar Foam::eConst::HaCell(const label celli) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return sMod_[Cp][celli]*T_->operator[](celli) + Hf_.value();
    }
    return
        sMod_[HDeparture][celli]
      + sMod_[Cp][celli]*T_->operator[](celli)
      + Hf_.value();
}


Foam::tmp<Foam::scalarField> Foam::eConst::HaPatch
(
    const label patchi
) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return
            sMod_[Cp].boundaryField()[patchi]
           *T_->boundaryField()[patchi]
          + Hf_.value();
    }
    return
        sMod_[HDeparture].boundaryField()[patchi]
      + sMod_[Cp].boundaryField()[patchi]
       *T_->boundaryField()[patchi]
      + Hf_.value();
}


Foam::tmp<Foam::scalarField> Foam::eConst::HaInternal() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return
            sMod_[Cp].primitiveField()
           *T_->primitiveField()
          + Hf_.value();
    }
    return
        sMod_[HDeparture].primitiveField()
      + sMod_[Cp].primitiveField()
       *T_->primitiveField()
      + Hf_.value();
}


Foam::tmp<Foam::volScalarField> Foam::eConst::HaGeometric() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return sMod_[Cp]()*T_->operator()() + Hf_.value();
    }
    return
        sMod_[Cp]()*T_->operator()() + Hf_.value() + sMod_[HDeparture]();
}


Foam::scalar Foam::eConst::HsCell(const label celli) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return sMod_[Cp][celli]*T_->operator[](celli);
    }
    return
        sMod_[HDeparture][celli]
      + sMod_[Cp][celli]*T_->operator[](celli);
}


Foam::tmp<Foam::scalarField> Foam::eConst::HsPatch
(
    const label patchi
) const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return
            sMod_[Cp].boundaryField()[patchi]
           *T_->boundaryField()[patchi];
    }
    return
        sMod_[HDeparture].boundaryField()[patchi]
      + sMod_[Cp].boundaryField()[patchi]*T_->boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::eConst::HsInternal() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return
            sMod_[Cp].primitiveField()*T_->primitiveField();
    }
    return
        sMod_[HDeparture].primitiveField()
      + sMod_[Cp].primitiveField()*T_->primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::eConst::HsGeometric() const
{
    if (sMod_(HDeparture) == nullptr)
    {
        return sMod_[Cp]()*T_->operator()();
    }
    return
        sMod_[Cp]()*T_->operator()() + sMod_[HDeparture]();
}


Foam::scalar Foam::eConst::SCell(const label celli) const
{
    const scalar S = sMod_[Cp][celli]*log(T_->operator[](celli)/Tstd);
    if (sMod_(SDeparture) == nullptr)
    {
        return S;
    }
    return S + sMod_[SDeparture][celli];
}


Foam::tmp<Foam::scalarField> Foam::eConst::SPatch
(
    const label patchi
) const
{
    tmp<scalarField> S
    (
        sMod_[Cp].boundaryField()[patchi]
       *log(T_->boundaryField()[patchi]/Tstd)
    );
    if (sMod_(SDeparture) == nullptr)
    {
        return S;
    }
    return S + sMod_[SDeparture].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::eConst::SInternal() const
{
    tmp<scalarField> S
    (
        sMod_[Cp].primitiveField()
       *log(T_->primitiveField()/Tstd)
    );
    if (sMod_(SDeparture) == nullptr)
    {
        return S;
    }
    return S + sMod_[SDeparture].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::eConst::SGeometric() const
{
    const dimensionedScalar TstdDim("Tstd", dimTemperature, Tstd);
    tmp<volScalarField> S(sMod_[Cp]()*log(T_->operator()()/TstdDim));
    if (sMod_(SDeparture) == nullptr)
    {
        return S;
    }
    return S + sMod_[SDeparture]();
}


bool Foam::eConst::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    Cv_.value() = dict_->lookup<scalar>("Cv");
    Hf_.value() = dict_->lookup<scalar>("Hf");

    return true;
}


// ************************************************************************* //
