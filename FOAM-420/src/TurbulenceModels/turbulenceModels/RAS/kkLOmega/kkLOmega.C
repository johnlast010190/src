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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2023 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "kkLOmega.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/bound/bound.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//defineTypeNameAndDebug(kkLOmega, 0);
//addToRunTimeSelectionTable(RASModel, kkLOmega, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::fv(const volScalarField& Ret) const
{
    return(1.0 - exp(-sqrt(Ret)/Av_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::fINT() const
{
    return
    (
        min
        (
            kt_/(Cint_*(kl_ + kt_ + this->kMin_)),
            dimensionedScalar("1.0", dimless, 1.0)
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::fSS(const volScalarField& Omega) const
{
    return(exp(-sqr(Css_*this->nu()*Omega/(kt_ + this->kMin_))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::Cmu(const volScalarField& S) const
{
    return(1.0/(A0_ + As_*(S/(omega_ + this->omegaMin_))));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::BetaTS(const volScalarField& ReOmega) const
{
    return(scalar(1) - exp(-sqr(max(ReOmega - CtsCrit_, scalar(0)))/Ats_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return
    (
        scalar(1)
      - exp
        (
            -CtauL_*ktL
          /
            (
                sqr
                (
                    lambdaEff*Omega
                  + dimensionedScalar
                    (
                        "ROOTVSMALL",
                        dimLength*inv(dimTime),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& ktS
) const
{
    return(fv*CmuStd_*sqrt(ktS)*lambdaEff);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::fOmega
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return
    (
        scalar(1)
      - exp
        (
           -0.41
           *pow4
            (
                lambdaEff
              / (
                    lambdaT
                  + dimensionedScalar
                    (
                        "ROTVSMALL",
                        lambdaT.dimensions(),
                        ROOTVSMALL
                    )
                )
            )
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::phiBP(const volScalarField& Omega) const
{
    return
    (
        min
        (
            max
            (
                kt_/this->nu()
             / (
                    Omega
                  + dimensionedScalar
                    (
                        "ROTVSMALL",
                        Omega.dimensions(),
                        ROOTVSMALL
                    )
                )
              - CbpCrit_,
                scalar(0)
            ),
            scalar(50.0)
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::phiNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    return
    (
        max
        (
            ReOmega
          - CnatCrit_
          / (
                fNatCrit + dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)
            ),
            scalar(0)
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kkLOmega<BasicTurbulenceModel>::D(const volScalarField& k) const
{
    return this->nu()*magSqr(fvc::grad(sqrt(k)));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kkLOmega<BasicTurbulenceModel>::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kkLOmega<BasicTurbulenceModel>::kkLOmega
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    As_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "As",
            this->coeffDict_,
            2.12
        )
    ),
    Av_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Av",
            this->coeffDict_,
            6.75
        )
    ),
    Abp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Abp",
            this->coeffDict_,
            0.6
        )
    ),
    Anat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anat",
            this->coeffDict_,
            200
        )
    ),
    Ats_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ats",
            this->coeffDict_,
            200
        )
    ),
    CbpCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CbpCrit",
            this->coeffDict_,
            1.2
        )
    ),
    Cnc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnc",
            this->coeffDict_,
            0.1
        )
    ),
    CnatCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatCrit",
            this->coeffDict_,
            1250
        )
    ),
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            this->coeffDict_,
            0.75
        )
    ),
    CtsCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit",
            this->coeffDict_,
            1000
        )
    ),
    CrNat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CrNat",
            this->coeffDict_,
            0.02
        )
    ),
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
            this->coeffDict_,
            3.4e-6
        )
    ),
    C12_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C12",
            this->coeffDict_,
            1.0e-10
        )
    ),
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
            this->coeffDict_,
            0.12
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
            this->coeffDict_,
            0.035
        )
    ),
    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            this->coeffDict_,
            1.5
        )
    ),
    CtauL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtauL",
            this->coeffDict_,
            4360
        )
    ),
    Cw1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw1",
            this->coeffDict_,
            0.44
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.92
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            0.3
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            this->coeffDict_,
            1.5
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
            this->coeffDict_,
            2.495
        )
    ),
    CmuStd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuStd",
            this->coeffDict_,
            0.09
        )
    ),
    Prtheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prtheta",
            this->coeffDict_,
            0.85
        )
    ),
    Sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmak",
            this->coeffDict_,
            1
        )
    ),
    Sigmaw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaw",
            this->coeffDict_,
            1.17
        )
    ),
    kt_
    (
        IOobject
        (
            IOobject::groupName("kt", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kl_
    (
        IOobject
        (
            IOobject::groupName("kl", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->db_
        ),
        kt_*omega_ + D(kl_) + D(kt_)
    )
{
    bound(kt_, this->kMin_);
    bound(kl_, this->kMin_);
    bound(omega_, this->omegaMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        // Evaluating nut_ is complex so start from the field read from file
        this->nut_.correctBoundaryConditions();

        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kkLOmega<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        A0_.readIfPresent(this->coeffDict());
        As_.readIfPresent(this->coeffDict());
        Av_.readIfPresent(this->coeffDict());
        Abp_.readIfPresent(this->coeffDict());
        Anat_.readIfPresent(this->coeffDict());
        Abp_.readIfPresent(this->coeffDict());
        Ats_.readIfPresent(this->coeffDict());
        CbpCrit_.readIfPresent(this->coeffDict());
        Cnc_.readIfPresent(this->coeffDict());
        CnatCrit_.readIfPresent(this->coeffDict());
        Cint_.readIfPresent(this->coeffDict());
        CtsCrit_.readIfPresent(this->coeffDict());
        CrNat_.readIfPresent(this->coeffDict());
        C11_.readIfPresent(this->coeffDict());
        C12_.readIfPresent(this->coeffDict());
        CR_.readIfPresent(this->coeffDict());
        CalphaTheta_.readIfPresent(this->coeffDict());
        Css_.readIfPresent(this->coeffDict());
        CtauL_.readIfPresent(this->coeffDict());
        Cw1_.readIfPresent(this->coeffDict());
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        CwR_.readIfPresent(this->coeffDict());
        Clambda_.readIfPresent(this->coeffDict());
        CmuStd_.readIfPresent(this->coeffDict());
        Prtheta_.readIfPresent(this->coeffDict());
        Sigmak_.readIfPresent(this->coeffDict());
        Sigmaw_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kkLOmega<BasicTurbulenceModel>::validate()
{}


template<class BasicTurbulenceModel>
void kkLOmega<BasicTurbulenceModel>::correct()
{
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;

    fv::options& fvOptions(fv::options::New(this->mesh_, omega_.db()));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField lambdaT(sqrt(kt_)/(omega_ + this->omegaMin_));

    const volScalarField& y(wallDist::New(this->mesh_).y());

    const volScalarField lambdaEff(min(Clambda_*y, lambdaT));

    const volScalarField fw
    (
        pow
        (
            lambdaEff
           /(lambdaT + dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2.0/3.0
        )
    );

    tmp<volTensorField> tgradU(fvc::grad(this->U_));
    const volTensorField& gradU = tgradU();

    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));

    const volScalarField S2(2.0*magSqr(dev(symm(gradU))));

    const volScalarField ktS(fSS(Omega)*fw*kt_);

    const volScalarField nuts
    (
        fv(sqr(fw)*kt_/this->nu()/(omega_ + this->omegaMin_))
       *fINT()
       *Cmu(sqrt(S2))*sqrt(ktS)*lambdaEff
    );
    const volScalarField Pkt(nuts*S2);

    const volScalarField ktL(kt_ - ktS);
    const volScalarField ReOmega(sqr(y)*Omega/this->nu());
    const volScalarField nutl
    (
        min
        (
            C11_*fTaul(lambdaEff, ktL, Omega)*Omega*sqr(lambdaEff)
           *sqrt(ktL)*lambdaEff/this->nu()
          + C12_*BetaTS(ReOmega)*ReOmega*sqr(y)*Omega
        ,
            0.5*(kl_ + ktL)/(sqrt(S2) + this->omegaMin_)
        )
    );

    const volScalarField Pkl(nutl*S2);

    const volScalarField alphaTEff
    (
        alphaT(lambdaEff, fv(sqr(fw)*kt_/this->nu()/(omega_ + this->omegaMin_)), ktS)
    );

    // By pass source term divided by kl_

    const dimensionedScalar fwMin("SMALL", dimless, ROOTVSMALL);

    const volScalarField Rbp
    (
        CR_*(1.0 - exp(-phiBP(Omega)()/Abp_))*omega_
       /(fw + fwMin)
    );

    const volScalarField fNatCrit(1.0 - exp(-Cnc_*sqrt(kl_)*y/this->nu()));

    // Natural source term divided by kl_
    const volScalarField Rnat
    (
        CrNat_*(1.0 - exp(-phiNAT(ReOmega, fNatCrit)/Anat_))*Omega
    );


    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(alphaTEff), omega_)
     ==
        Cw1_*Pkt*alpha*rho*omega_/(kt_ + this->kMin_)
      - fvm::SuSp
        (
            (1.0 - CwR_/(fw + fwMin))*kl_*(Rbp + Rnat)*alpha*rho/(kt_ + this->kMin_),
            omega_
        )
      - fvm::Sp(Cw2_*sqr(fw)*alpha*rho*omega_, omega_)
      + (
            Cw3_*fOmega(lambdaEff, lambdaT)*alpha*rho*alphaTEff*sqr(fw)*sqrt(kt_)
        )()()/pow3(y.internalField())
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());

    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    const volScalarField Dl(D(kl_));

    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(alpha, rho, kl_)
      + fvm::div(alphaRhoPhi, kl_)
      - fvm::laplacian(alpha*rho*this->nu(), kl_)
     ==
        Pkl*alpha*rho
      - fvm::Sp((Rbp + Rnat + Dl/(kl_ + this->kMin_))*alpha*rho, kl_)
      + fvOptions(alpha, rho, kl_)
    );

    klEqn.ref().relax();
    fvOptions.constrain(klEqn.ref());

    solve(klEqn);
    fvOptions.correct(kl_);
    bound(kl_, this->kMin_);


    const volScalarField Dt(D(kt_));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ktEqn
    (
        fvm::ddt(alpha, rho, kt_)
      + fvm::div(alphaRhoPhi, kt_)
      - fvm::laplacian(alpha*rho*DkEff(alphaTEff), kt_)
     ==
        Pkt*alpha*rho
      + (Rbp + Rnat)*alpha*rho*kl_
      - fvm::Sp((omega_ + Dt/(kt_+ this->kMin_))*alpha*rho, kt_)
      + fvOptions(alpha, rho, kt_)
    );

    ktEqn.ref().relax();
    fvOptions.constrain(ktEqn.ref());

    solve(ktEqn);
    fvOptions.correct(kt_);
    bound(kt_, this->kMin_);


    // Update total fluctuation kinetic energy dissipation rate
    epsilon_ = kt_*omega_ + Dl + Dt;
    bound(epsilon_, this->epsilonMin_);


    // Re-calculate turbulent viscosity
    this->nut_ = nuts + nutl;
    this->nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
