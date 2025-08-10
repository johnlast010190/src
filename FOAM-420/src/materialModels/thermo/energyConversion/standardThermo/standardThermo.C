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

#include "standardThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(standardThermo, 0);
    addToRunTimeSelectionTable(materialModel, standardThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::standardThermo::standardThermo
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
    dep_.setSize(6);
    standardThermo::read();
}


Foam::autoPtr<Foam::standardThermo> Foam::standardThermo::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<standardThermo>
    (
        new standardThermo
        (
            obr,
            dict,
            phaseName,
            specieName,
            name
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::standardThermo::~standardThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::standardThermo::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    // Model names
    const word& CvName = CvModel::typeName;
    const word& CpName = CpModel::typeName;
    const word& CpMCvName = CpMCvModel::typeName;
    const word& gammaName = gammaModel::typeName;
    const word& EsName = EsModel::typeName;
    const word& EaName = EaModel::typeName;
    const word& rhoName = rhoModel::typeName;
    const word& HsName = HsModel::typeName;
    const word& HaName = HaModel::typeName;
    const word& GName = GModel::typeName;
    const word& AName = AModel::typeName;
    const word& SName = SModel::typeName;

    if (modelName == CvName)
    {
        sMod_.set(Cp, models[CpName]);
        sMod_.set(CpMCv, models[CpMCvName]);
        dep_[0].model = models[CvName];
        dep_[0].dependencies.setSize(2);
        dep_[0].dependencies.set(0, models[CpName]);
        dep_[0].dependencies.set(1, models[CpMCvName]);
    }
    else if (modelName == gammaName)
    {
        sMod_.set(Cp, models[CpName]);
        sMod_.set(CpMCv, models[CpMCvName]);
        dep_[1].model = models[gammaName];
        dep_[1].dependencies.setSize(2);
        dep_[1].dependencies.set(0, models[CpName]);
        dep_[1].dependencies.set(1, models[CpMCvName]);
    }
    else if (modelName == EsName)
    {
        sMod_.set(rho, models[rhoName]);
        sMod_.set(Hs, models[HsName]);
        dep_[2].model = models[EsName];
        dep_[2].dependencies.setSize(2);
        dep_[2].dependencies.set(0, models[rhoName]);
        dep_[2].dependencies.set(1, models[HsName]);
    }
    else if (modelName == EaName)
    {
        sMod_.set(rho, models[rhoName]);
        sMod_.set(Ha, models[HaName]);
        dep_[3].model = models[EaName];
        dep_[3].dependencies.setSize(2);
        dep_[3].dependencies.set(0, models[rhoName]);
        dep_[3].dependencies.set(1, models[HaName]);
    }
    else if (modelName == GName)
    {
        sMod_.set(Ha, models[HaName]);
        sMod_.set(S, models[SName]);
        dep_[4].model = models[GName];
        dep_[4].dependencies.setSize(2);
        dep_[4].dependencies.set(0, models[HaName]);
        dep_[4].dependencies.set(1, models[SName]);
    }
    else if (modelName == AName)
    {
        sMod_.set(Ea, models[EaName]);
        sMod_.set(S, models[SName]);
        dep_[5].model = models[AName];
        dep_[5].dependencies.setSize(2);
        dep_[5].dependencies.set(0, models[EaName]);
        dep_[5].dependencies.set(1, models[SName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::standardThermo::castScalarModel
(
    const word& modelName
)
{
    if (modelName == AModel::typeName)
    {
        return dynamic_cast<AModel*>(this);
    }
    else if (modelName == CvModel::typeName)
    {
        return dynamic_cast<CvModel*>(this);
    }
    else if (modelName == EaModel::typeName)
    {
        return dynamic_cast<EaModel*>(this);
    }
    else if (modelName == EsModel::typeName)
    {
        return dynamic_cast<EsModel*>(this);
    }
    else if (modelName == gammaModel::typeName)
    {
        return dynamic_cast<gammaModel*>(this);
    }
    else if (modelName == GModel::typeName)
    {
        return dynamic_cast<GModel*>(this);
    }
    return nullptr;
}


Foam::scalar Foam::standardThermo::CvCell(const label celli) const
{
    return sMod_[Cp][celli] - sMod_[CpMCv][celli];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::CvPatch
(
    const label patchi
) const
{
    return
        sMod_[Cp].boundaryField()[patchi]
      - sMod_[CpMCv].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::CvInternal() const
{
    return sMod_[Cp].primitiveField() - sMod_[CpMCv].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::standardThermo::CvGeometric() const
{
    return sMod_[Cp]() - sMod_[CpMCv]();
}


Foam::scalar Foam::standardThermo::gammaCell(const label celli) const
{
    const scalar CpNew = sMod_[Cp][celli];
    return CpNew/(CpNew - sMod_[CpMCv][celli]);
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::gammaPatch
(
    const label patchi
) const
{
    const tmp<scalarField> tCp(sMod_[Cp].boundaryField()[patchi]);
    return tCp.ref()/(tCp.ref() - sMod_[CpMCv].boundaryField()[patchi]);
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::gammaInternal() const
{
    const tmp<scalarField> tCp(sMod_[Cp].primitiveField());
    return tCp.ref()/(tCp.ref() - sMod_[CpMCv].primitiveField());
}


Foam::tmp<Foam::volScalarField> Foam::standardThermo::gammaGeometric() const
{
    const tmp<volScalarField> tCp(sMod_[Cp]());
    return tCp.ref()/(tCp.ref() - sMod_[CpMCv]());
}


Foam::scalar Foam::standardThermo::EsCell(const label celli) const
{
    return
        sMod_[Hs][celli] - p_->operator[](celli)/sMod_[rho][celli];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::EsPatch
(
    const label patchi
) const
{
    return
        sMod_[Hs].boundaryField()[patchi]
      - p_->boundaryField()[patchi]
       /sMod_[rho].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::EsInternal() const
{
    return
          sMod_[Hs].primitiveField()
        - p_->primitiveField()
         /sMod_[rho].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::standardThermo::EsGeometric() const
{
    return sMod_[Hs]() - p_->operator()()/sMod_[rho]();
}


Foam::scalar Foam::standardThermo::EaCell(const label celli) const
{
    return
        sMod_[Ha][celli] - p_->operator[](celli)/sMod_[rho][celli];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::EaPatch
(
    const label patchi
) const
{
    return
        sMod_[Ha].boundaryField()[patchi]
      - p_->boundaryField()[patchi]
       /sMod_[rho].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::EaInternal() const
{
    return
        sMod_[Ha].primitiveField()
      - p_->primitiveField()
       /sMod_[rho].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::standardThermo::EaGeometric() const
{
    return sMod_[Ha]() - p_->operator()()/sMod_[rho]();
}


Foam::scalar Foam::standardThermo::GCell(const label celli) const
{
    return
        sMod_[Ha][celli] - T_->operator[](celli)*sMod_[S][celli];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::GPatch
(
    const label patchi
) const
{
    return
        sMod_[Ha].boundaryField()[patchi]
      - T_->boundaryField()[patchi]
       *sMod_[S].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::GInternal() const
{
    return
        sMod_[Ha].primitiveField()
      - T_->primitiveField()
       *sMod_[S].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::standardThermo::GGeometric() const
{
    return sMod_[Ha]() - T_->operator()()*sMod_[S]();
}


Foam::scalar Foam::standardThermo::ACell(const label celli) const
{
    return sMod_[Ea][celli] - T_->operator[](celli)*sMod_[S][celli];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::APatch
(
    const label patchi
) const
{
    return
        sMod_[Ea].boundaryField()[patchi]
      - T_->boundaryField()[patchi]
       *sMod_[S].boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::standardThermo::AInternal() const
{
    return
        sMod_[Ea].primitiveField()
      - T_->primitiveField()
       *sMod_[S].primitiveField();
}


Foam::tmp<Foam::volScalarField> Foam::standardThermo::AGeometric() const
{
    return sMod_[Ea]() - T_->operator()()*sMod_[S]();
}


bool Foam::standardThermo::read()
{
    p_ = constructOrReturnRefFieldPtr<scalar>("p");
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    return true;
}


// ************************************************************************* //
