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
    (c) 2014-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "PengRobinsonGasLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PengRobinsonGasLaw, 0);
    addToRunTimeSelectionTable(materialModel, PengRobinsonGasLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PengRobinsonGasLaw::PengRobinsonGasLaw
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
    PengRobinsonGasLaw::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::PengRobinsonGasLaw::incompressible() const
{
    if (p_->isConst())
    {
        return true;
    }
    return false;
}


bool Foam::PengRobinsonGasLaw::isochoric() const
{
    if (T_->isConst() && p_->isConst())
    {
        return true;
    }
    return false;
}


void Foam::PengRobinsonGasLaw::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    const word tName = materialTables_.tableName(phaseName_, specieName_);

    // Model names
    const word& RName = RModel::typeName;
    const word& WName = WModel::typeName;
    const word& SDepartureName = SDepartureModel::typeName;
    const word& HDepartureName = HDepartureModel::typeName;
    const word& CpDepartureName = CpDepartureModel::typeName;
    const word& psiName = psiModel::typeName;
    const word& rhoName = rhoModel::typeName;
    const word& CpMCvName = CpMCvModel::typeName;
    const word& ZName = ZModel::typeName;

    if (modelName == rhoName)
    {
        sMod_.set(Rind, models[RName]);
        sMod_.set(Zind, models[ZName]);

        // Dependencies
        dep_[0].model = models[rhoName];
        dep_[0].dependencies.setSize(2);
        dep_[0].dependencies.set(0, models[RName]);
        dep_[0].dependencies.set(1, models[ZName]);
    }
    else if (modelName == HDepartureName)
    {
        sMod_.set(Rind, models[RName]);
        sMod_.set(Zind, models[ZName]);

        // Dependencies
        dep_[1].model = models[HDepartureName];
        dep_[1].dependencies.setSize(2);
        dep_[1].dependencies.set(0, models[RName]);
        dep_[1].dependencies.set(1, models[ZName]);
    }
    else if (modelName == CpDepartureName)
    {
        sMod_.set(Rind, models[RName]);
        sMod_.set(Wind, models[WName]);

        // Dependencies
        dep_[2].model = models[CpDepartureName];
        dep_[2].dependencies.setSize(2);
        dep_[2].dependencies.set(0, models[RName]);
        dep_[2].dependencies.set(1, models[WName]);
    }
    else if (modelName == SDepartureName)
    {
        sMod_.set(Rind, models[RName]);
        sMod_.set(Zind, models[ZName]);

        // Dependencies
        dep_[3].model = models[SDepartureName];
        dep_[3].dependencies.setSize(2);
        dep_[3].dependencies.set(0, models[RName]);
        dep_[3].dependencies.set(1, models[ZName]);
    }
    else if (modelName == psiName)
    {
        sMod_.set(Rind, models[RName]);
        sMod_.set(Zind, models[ZName]);

        // Dependencies
        dep_[4].model = models[psiName];
        dep_[4].dependencies.setSize(2);
        dep_[4].dependencies.set(0, models[RName]);
        dep_[4].dependencies.set(1, models[ZName]);
    }
    else if (modelName == CpMCvName)
    {
        sMod_.set(Rind, models[RName]);
        sMod_.set(Zind, models[ZName]);

        // Dependencies
        dep_[5].model = models[CpMCvName];
        dep_[5].dependencies.setSize(2);
        dep_[5].dependencies.set(0, models[RName]);
        dep_[5].dependencies.set(1, models[ZName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::PengRobinsonGasLaw::castScalarModel
(
    const word& modelName
)
{
    if (modelName == rhoModel::typeName)
    {
        return dynamic_cast<rhoModel*>(this);
    }
    else if (modelName == HDepartureModel::typeName)
    {
        return dynamic_cast<HDepartureModel*>(this);
    }
    else if (modelName == CpDepartureModel::typeName)
    {
        return dynamic_cast<CpDepartureModel*>(this);
    }
    else if (modelName == SDepartureModel::typeName)
    {
        return dynamic_cast<SDepartureModel*>(this);
    }
    else if (modelName == psiModel::typeName)
    {
        return dynamic_cast<psiModel*>(this);
    }
    else if (modelName == ZModel::typeName)
    {
        return dynamic_cast<ZModel*>(this);
    }
    else if (modelName == CpMCvModel::typeName)
    {
        return dynamic_cast<CpMCvModel*>(this);
    }
    return nullptr;
}


Foam::scalar Foam::PengRobinsonGasLaw::rhoCell(const label celli) const
{
    return
        p_->operator[](celli)
        /(sMod_[Zind][celli]*sMod_[Rind][celli]*T_->operator[](celli));
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::rhoPatch
(
    const label patchi
) const
{
    return
        p_->boundaryField()[patchi]
       /(
            sMod_[Zind].boundaryField()[patchi]
           *sMod_[Rind].boundaryField()[patchi]
           *T_->boundaryField()[patchi]
        );
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::rhoInternal() const
{
    return
        p_->primitiveField()
       /(
            sMod_[Zind].primitiveField()
           *sMod_[Rind].primitiveField()
           *T_->primitiveField()
        );
}


Foam::tmp<Foam::volScalarField>
Foam::PengRobinsonGasLaw::rhoGeometric() const
{
    return p_->operator()()/(sMod_[Zind]()*sMod_[Rind]()*T_->operator()());
}


Foam::scalar Foam::PengRobinsonGasLaw::HDepartureCell(const label celli) const
{
    const scalar Pr = p_->operator[](celli)/Pc_;
    const scalar Tr = T_->operator[](celli)/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar alpha = sqr(1 + kappa_*(1 - sqrt(Tr)));
    const scalar Z = sMod_[Zind][celli];

    return
        sMod_[Rind][celli]
       *Tc_
       *(
           Tr*(Z - 1)
         - 2.078*(1 + kappa_)*sqrt(alpha)
          *log((Z + 2.414*B)/(Z - 0.414*B))
        );
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::HDeparturePatch
(
    const label patchi
) const
{
    tmp<scalarField> Pr(p_->boundaryField()[patchi]/Pc_);
    tmp<scalarField> Tr(T_->boundaryField()[patchi]/Tc_);
    tmp<scalarField> B(0.07780*Pr/Tr.ref());
    tmp<scalarField> alpha(sqr(1 + kappa_*(1 - sqrt(Tr.ref()))));
    tmp<scalarField> Z(sMod_[Zind].boundaryField()[patchi]);

    return
        sMod_[Rind].boundaryField()[patchi]
       *Tc_
       *(
           Tr.ref()*(Z.ref() - 1)
         - 2.078*(1 + kappa_)*sqrt(alpha)
          *log((Z.ref() + 2.414*B.ref())/(Z.ref() - 0.414*B.ref()))
        );
}


Foam::tmp<Foam::scalarField>
Foam::PengRobinsonGasLaw::HDepartureInternal() const
{
    tmp<scalarField> Pr(p_->primitiveField()/Pc_);
    tmp<scalarField> Tr(T_->primitiveField()/Tc_);
    tmp<scalarField> B(0.07780*Pr/Tr.ref());
    tmp<scalarField> alpha(sqr(1 + kappa_*(1 - sqrt(Tr.ref()))));
    tmp<scalarField> Z(sMod_[Zind].primitiveField());

    return
        sMod_[Rind].primitiveField()
       *Tc_
       *(
           Tr.ref()*(Z.ref() - 1)
         - 2.078*(1 + kappa_)*sqrt(alpha)
          *log((Z.ref() + 2.414*B.ref())/(Z.ref() - 0.414*B.ref()))
        );
}


Foam::tmp<Foam::volScalarField>
Foam::PengRobinsonGasLaw::HDepartureGeometric() const
{
    tmp<volScalarField> tHDeparture
    (
        new volScalarField
        (
            IOobject
            (
                "HDeparture",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar
            (
                "HDeparture",
                HDepartureModel::modelDims,
                Zero
            )
        )
    );
    volScalarField& HDeparture = tHDeparture.ref();

    HDeparture.primitiveFieldRef() = HDepartureInternal();
    forAll(HDeparture.boundaryField(), patchi)
    {
        HDeparture.boundaryFieldRef()[patchi].forceAssign(
            HDeparturePatch(patchi)
        );
    }

    return tHDeparture;
}


Foam::scalar Foam::PengRobinsonGasLaw::CpDepartureCell(const label celli) const
{
    const scalar Tr = T_->operator[](celli)/Tc_;
    const scalar alpha = sqr(1 + kappa_*(1 - sqrt(Tr)));
    const scalar A =
        a_*alpha*p_->operator[](celli)
       /sqr(RR*T_->operator[](celli));
    const scalar B =
        b_*p_->operator[](celli)/(RR*T_->operator[](celli));
    const scalar Z = sMod_[Zind][celli];
    const scalar ap =
        kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T_->operator[](celli)*Tc_));
    const scalar app =
        kappa_*a_*(1 + kappa_)/(2*sqrt(pow3(T_->operator[](celli))*Tc_));
    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    const scalar root2 = sqrt(2.0);

    return
    (
        app
       *(T_->operator[](celli)/(2*root2*b_))
       *log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
      + RR*sqr(M - N)/(sqr(M) - 2*A*(Z + B))
      - RR
    )/sMod_[Wind][celli];
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::CpDeparturePatch
(
    const label patchi
) const
{
    const tmp<scalarField> Tr(T_->boundaryField()[patchi]/Tc_);
    const tmp<scalarField> alpha(sqr(1 + kappa_*(1 - sqrt(Tr))));
    const tmp<scalarField> A
    (
        a_*alpha*p_->boundaryField()[patchi]
       /sqr(RR*T_->boundaryField()[patchi])
    );
    const tmp<scalarField> B
    (
        b_*p_->boundaryField()[patchi]
        /(RR*T_->boundaryField()[patchi])
    );
    const tmp<scalarField> Z(sMod_[Zind].boundaryField()[patchi]);
    const tmp<scalarField> ap
    (
        kappa_*a_*(kappa_/Tc_ - (1 + kappa_)
       /sqrt(T_->boundaryField()[patchi]*Tc_))
    );
    const tmp<scalarField> app
    (
        kappa_*a_*(1 + kappa_)/(2*sqrt(pow3(T_->boundaryField()[patchi])*Tc_))
    );
    const tmp<scalarField> M
    (
        (sqr(Z.ref()) + 2*B.ref()*Z.ref() - sqr(B.ref()))/(Z.ref() - B.ref())
    );
    const tmp<scalarField> N(ap*B.ref()/(b_*RR));

    const scalar root2 = sqrt(2.0);

    return
    (
        app
       *(T_->boundaryField()[patchi]/(2*root2*b_))
       *log((Z.ref() + (root2 + 1)*B.ref())/(Z.ref() - (root2 - 1)*B.ref()))
      + RR*sqr(M.ref() - N)/(sqr(M.ref()) - 2*A*(Z.ref() + B.ref()))
      - RR
    )/sMod_[Wind][0];
}


Foam::tmp<Foam::scalarField>
Foam::PengRobinsonGasLaw::CpDepartureInternal() const
{
    const tmp<scalarField> Tr(T_->primitiveField()/Tc_);
    const tmp<scalarField> alpha(sqr(1 + kappa_*(1 - sqrt(Tr))));
    const tmp<scalarField> A
    (
        a_*alpha*p_->primitiveField()/sqr(RR*T_->primitiveField())
    );
    const tmp<scalarField> B
    (
        b_*p_->primitiveField()/(RR*T_->primitiveField())
    );
    const tmp<scalarField> Z(sMod_[Zind].primitiveField());
    const tmp<scalarField> ap
    (
        kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T_->primitiveField()*Tc_))
    );
    const tmp<scalarField> app
    (
        kappa_*a_*(1 + kappa_)/(2*sqrt(pow3(T_->primitiveField())*Tc_))
    );
    const tmp<scalarField> M
    (
        (sqr(Z.ref()) + 2*B.ref()*Z.ref() - sqr(B.ref()))/(Z.ref() - B.ref())
    );
    const tmp<scalarField> N(ap*B.ref()/(b_*RR));

    const scalar root2 = sqrt(2.0);

    return
    (
        app
       *(T_->primitiveField()/(2*root2*b_))
       *log((Z.ref() + (root2 + 1)*B.ref())/(Z.ref() - (root2 - 1)*B.ref()))
      + RR*sqr(M.ref() - N)/(sqr(M.ref()) - 2*A*(Z.ref() + B.ref()))
      - RR
    )/sMod_[Wind][0];
}


Foam::tmp<Foam::volScalarField>
Foam::PengRobinsonGasLaw::CpDepartureGeometric() const
{
    tmp<volScalarField> tCpDeparture
    (
        new volScalarField
        (
            IOobject
            (
                "CpDeparture",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar
            (
                "CpDeparture",
                CpDepartureModel::modelDims,
                Zero
            )
        )
    );
    volScalarField& CpDeparture = tCpDeparture.ref();

    CpDeparture.primitiveFieldRef() = CpDepartureInternal();
    forAll(CpDeparture.boundaryField(), patchi)
    {
        CpDeparture.boundaryFieldRef()[patchi].forceAssign(
            CpDeparturePatch(patchi)
        );
    }

    return tCpDeparture;
}


Foam::scalar Foam::PengRobinsonGasLaw::SDepartureCell(const label celli) const
{
    const scalar Pr = p_->operator[](celli)/Pc_;
    const scalar Tr = T_->operator[](celli)/Tc_;
    const scalar B = 0.07780*Pr/Tr;
    const scalar Z = sMod_[Zind][celli];

    return
        sMod_[Rind][celli]
       *(
          - log(p_->operator[](celli)/Pstd)
          + (
                log(Z - B)
              - 2.078*kappa_*((1 + kappa_)/sqrt(Tr) - kappa_)
               *log((Z + 2.414*B)/(Z - 0.414*B))
            )
        );
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::SDeparturePatch
(
    const label patchi
) const
{
    const tmp<scalarField> Pr(p_->boundaryField()[patchi]/Pc_);
    const tmp<scalarField> Tr(T_->boundaryField()[patchi]/Tc_);
    const tmp<scalarField> B(0.07780*Pr/Tr.ref());
    const tmp<scalarField> Z(sMod_[Zind].boundaryField()[patchi]);

    return
        sMod_[Rind][0]
       *(
          - log(p_->boundaryField()[patchi]/Pstd)
          + (
                log(Z.ref() - B.ref())
              - 2.078*kappa_*((1 + kappa_)/sqrt(Tr.ref()) - kappa_)
               *log((Z.ref() + 2.414*B.ref())/(Z.ref() - 0.414*B.ref()))
            )
        );
}


Foam::tmp<Foam::scalarField>
Foam::PengRobinsonGasLaw::SDepartureInternal() const
{
    const tmp<scalarField> Pr(p_->primitiveField()/Pc_);
    const tmp<scalarField> Tr(T_->primitiveField()/Tc_);
    const tmp<scalarField> B(0.07780*Pr/Tr.ref());
    const tmp<scalarField> Z(sMod_[Zind].primitiveField());

    return
        sMod_[Rind][0]
       *(
          - log(p_->primitiveField()/Pstd)
          + (
                log(Z.ref() - B.ref())
              - 2.078*kappa_*((1 + kappa_)/sqrt(Tr.ref()) - kappa_)
               *log((Z.ref() + 2.414*B.ref())/(Z.ref() - 0.414*B.ref()))
            )
        );
}


Foam::tmp<Foam::volScalarField>
Foam::PengRobinsonGasLaw::SDepartureGeometric() const
{
    tmp<volScalarField> tSDeparture
    (
        new volScalarField
        (
            IOobject
            (
                "SDeparture",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar
            (
                "SDeparture",
                SDepartureModel::modelDims,
                Zero
            )
        )
    );
    volScalarField& SDeparture = tSDeparture.ref();

    SDeparture.primitiveFieldRef() = SDepartureInternal();
    forAll(SDeparture.boundaryField(), patchi)
    {
        SDeparture.boundaryFieldRef()[patchi].forceAssign(
            SDeparturePatch(patchi)
        );
    }

    return tSDeparture;
}


Foam::scalar Foam::PengRobinsonGasLaw::psiCell(const label celli) const
{
    return 1.0/(sMod_[Zind][celli]*sMod_[Rind][0]*T_->operator[](celli));
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::psiPatch
(
    const label patchi
) const
{
    return
        1.0
       /(
           sMod_[Zind].boundaryField()[patchi]
          *sMod_[Rind][0]*T_->boundaryField()[patchi]
        );
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::psiInternal() const
{
    return
        1.0
       /(
           sMod_[Zind].primitiveField()
          *sMod_[Rind][0]*T_->primitiveField()
        );
}


Foam::tmp<Foam::volScalarField>
Foam::PengRobinsonGasLaw::psiGeometric() const
{

    return 1.0/(sMod_[Zind]()*sMod_[Rind]()*T_->operator()());
}


Foam::scalar Foam::PengRobinsonGasLaw::CpMCvCell(const label celli) const
{
    const scalar Tr = T_->operator[](celli)/Tc_;
    const scalar alpha = sqr(1 + kappa_*(1 - sqrt(Tr)));
    const scalar A = alpha*a_*p_->operator[](celli)/sqr(RR*T_->operator[](celli));
    const scalar B = b_*p_->operator[](celli)/(RR*T_->operator[](celli));
    const scalar Z = sMod_[Zind][celli];
    const scalar ap =
        kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T_->operator[](celli)*Tc_));
    const scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
    const scalar N = ap*B/(b_*RR);

    return sMod_[Rind][0]*sqr(M - N)/(sqr(M) - 2*A*(Z + B));
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::CpMCvPatch
(
    const label patchi
) const
{
    const tmp<scalarField> Tr(T_->boundaryField()[patchi]/Tc_);
    const tmp<scalarField> alpha = sqr(1 + kappa_*(1 - sqrt(Tr.ref())));

    const tmp<scalarField> A
    (
        alpha*a_*p_->boundaryField()[patchi]
       /sqr(RR*T_->boundaryField()[patchi])
    );
    const tmp<scalarField> B
    (
        b_*p_->boundaryField()[patchi]
        /(RR*T_->boundaryField()[patchi])
    );
    const tmp<scalarField> Z(sMod_[Zind].boundaryField()[patchi]);
    const tmp<scalarField> ap
    (
        kappa_*a_*(kappa_/Tc_ - (1 + kappa_)
       /sqrt(T_->boundaryField()[patchi]*Tc_))
    );
    const tmp<scalarField> M
    (
        (sqr(Z.ref()) + 2*B.ref()*Z.ref() - sqr(B.ref()))/(Z.ref() - B.ref())
    );
    const tmp<scalarField> N(ap*B.ref()/(b_*RR));

    return
        sMod_[Rind][0]*sqr(M.ref() - N.ref())
       /(sqr(M.ref()) - 2*A.ref()*(Z.ref() + B.ref()));
}


Foam::tmp<Foam::scalarField> Foam::PengRobinsonGasLaw::CpMCvInternal() const
{
    const tmp<scalarField> Tr(T_->primitiveField()/Tc_);
    const tmp<scalarField> alpha = sqr(1 + kappa_*(1 - sqrt(Tr.ref())));
    const tmp<scalarField> A
    (
        alpha*a_*p_->primitiveField()/sqr(RR*T_->primitiveField())
    );
    const tmp<scalarField> B
    (
        b_*p_->primitiveField()/(RR*T_->primitiveField())
    );
    const tmp<scalarField> Z(sMod_[Zind].primitiveField());
    const tmp<scalarField> ap
    (
        kappa_*a_*(kappa_/Tc_ - (1 + kappa_)/sqrt(T_->primitiveField()*Tc_))
    );
    const tmp<scalarField> M
    (
        (sqr(Z.ref()) + 2*B.ref()*Z.ref() - sqr(B.ref()))/(Z.ref() - B.ref())
    );
    const tmp<scalarField> N(ap*B.ref()/(b_*RR));

    return
        sMod_[Rind][0]*sqr(M.ref() - N.ref())
       /(sqr(M.ref()) - 2*A.ref()*(Z.ref() + B.ref()));
}


Foam::tmp<Foam::volScalarField>
Foam::PengRobinsonGasLaw::CpMCvGeometric() const
{
    tmp<volScalarField> tCpMCv
    (
        new volScalarField
        (
            IOobject
            (
                "CpMCv",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar
            (
                "CpMCv",
                CpMCvModel::modelDims,
                Zero
            )
        )
    );
    volScalarField& CpMCv = tCpMCv.ref();

    CpMCv.primitiveFieldRef() = CpMCvInternal();
    forAll(CpMCv.boundaryField(), patchi)
    {
        CpMCv.boundaryFieldRef()[patchi].forceAssign(CpMCvPatch(patchi));
    }

    return tCpMCv;
}


bool Foam::PengRobinsonGasLaw::read()
{
    p_ = constructOrReturnRefFieldPtr<scalar>("p");
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    if (sMod_(Rind) != nullptr)
    {
        sMod_[Rind].read();
        R_ = sMod_[Rind][0];
    }
    Tc_= dict_->lookup<scalar>("Tc");
    Vc_ = dict_->lookup<scalar>("Vc");
    Pc_ = dict_->lookup<scalar>("Pc");
    omega_ = dict_->lookup<scalar>("omega");
    Zc_ = Pc_*Vc_/(RR*Tc_);
    a_ = 0.45724*sqr(RR*Tc_)/Pc_;
    b_ = 0.07780*RR*Tc_/Pc_;
    kappa_ = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);

    return true;
}


// ************************************************************************* //
