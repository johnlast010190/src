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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "ellipticBlendingLagKE.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/bound/bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> ellipticBlendingLagKE<BasicTurbulenceModel>::Ts() const
{
    //  ( Teddy + Tkol )
    return sqrt(k_*k_/sqr(epsilon_ + dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL)) +
                CT_*CT_*(max(this->nu(), dimensionedScalar("zero", this->nu()().dimensions(), 0.0))
               /(epsilon_ + dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL))));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> ellipticBlendingLagKE<BasicTurbulenceModel>::TsLim
(
    const volTensorField& gradU
) const
{
    //  ( Treal )
    return CtLim_/(sqrt(3*2*magSqr(dev(symm(gradU))))*Cmu_*v2k_ +
               dimensionedScalar("TinvMin", dimensionSet(0,0,-1,0,0,0,0), SMALL));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> ellipticBlendingLagKE<BasicTurbulenceModel>::Ls() const
{
    return CL_*sqrt((pow(k_,3.0)/
           (pow(epsilon_+ dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL), 2.0)))
          + Ceta_*Ceta_*sqrt(pow(this->nu(),3.0)/(epsilon_ +
               dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL))));
}

template<class BasicTurbulenceModel>
void ellipticBlendingLagKE<BasicTurbulenceModel>::correctNut
(
    const volTensorField& gradU
)
{
    const volScalarField Ts(this->Ts());
    const volScalarField TsLim(this->TsLim(gradU));

    this->nut_ = this->Cmu_*k_*v2k_*min(Ts, TsLim);
    nutOverk_ = this->Cmu_*v2k_*min(Ts, TsLim);

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_, this->nut_.db()).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void ellipticBlendingLagKE<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2
)
{
    correctNut(fvc::grad(this->U_));
}

template<class BasicTurbulenceModel>
void ellipticBlendingLagKE<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}

template<class BasicTurbulenceModel>
void ellipticBlendingLagKE<BasicTurbulenceModel>::solveTurb
(
    const volTensorField& gradU
)
{
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_, U.db()));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    const volSymmTensorField S("Strain", dev(symm(gradU)));
    volScalarField S2(2*magSqr(symm(gradU)));
    const volScalarField sqrtS2(sqrt(S2));
    const volScalarField G(this->GName(), nut*S2);

    const volScalarField Goverk("GoverK", nutOverk_*S2);

    const volScalarField L2(type() + ":L2", sqr(Ls()));

    const volSymmTensorField Sbar ("Sbar", S/(sqrtS2+dimensionedScalar("small", dimensionSet(0,0,-1,0,0), SMALL)));
    const volSymmTensorField dSijDt (fvc::DDt(this->phi(), Sbar));


    volTensorField Wshur
    (
        "Wshur",
        (dSijDt & Sbar) - (Sbar & dSijDt)
    );

    surfaceTensorField WshurS( fvc::interpolate(Wshur) );
    Wshur = fvc::average(WshurS);

    const volTensorField W(skew(gradU) - Wshur);
    const volTensorField SW(S + W);

//    const volScalarField Ceps2star
//    (
//        "Ceps2star",
//        Ceps2_ + pow(a_,3)*(Ceps4_ - Ceps2_)*tanh(pow(mag(fvc::div(nut/sigmaK_*fvc::grad(k_))/epsilon_),1.5))
//    );
    const volVectorField& n_(wallDist::New(this->mesh_).n());
//    const volVectorField grada(fvc::grad(a_));
//    const volVectorField n_(grada/(mag(grada) + dimensionedScalar("nsmall", grada.dimensions(), 1.)));
    const volScalarField E
    (
        "E", 2.0*this->nu()*nut*sqr(fvc::div(mag(2*S & n_ ) * n_, "div(Sn*n)" ))
    );

    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1_*alpha*rho*Goverk*epsilon_
      - fvm::SuSp(Ceps2_*alpha*rho*epsilon_/k_, epsilon_)
//      - fvm::SuSp(Ceps2star*alpha*rho*epsilon_/k_, epsilon_)  // for variable Ceps2
      + Ceps3_*pow(1.0 - a_,3.0)*E*alpha*rho
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    volScalarField EpsCoeff
    (
        "EpsCoeff",
        epsilon_/k_
    );
    scalar unityCoef = 1;

    const volScalarField Sk
    (
        "Sk",
        unityCoef*G
    );

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*Sk
      - fvm::SuSp(2.0/3.0*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*EpsCoeff, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    tmp<fvScalarMatrix> aEqn
    (
       -fvm::laplacian(a_)
        ==
        1.0/L2
       -fvm::Sp(1.0/L2, a_)
    );

    aEqn.ref().relax();
    fvOptions.constrain(aEqn.ref());
    solve(aEqn);
    fvOptions.correct(a_);
    bound(a_, aMin_);

    // Turbulence stress normal to streamlines equation
    correctNut(gradU);
    const volScalarField Ts(this->Ts());
    const volTensorField Omikron(2*W/(sqrt(SW && SW)+dimensionedScalar("small", dimensionSet(0,0,-1,0,0), SMALL)));
    const volTensorField b
    (
        -2.0*nutOverk_*(S - 0.888*(Omikron&S))
    );

    const volScalarField fmu
    (
       "fmu",
       (sqrtS2*k_/epsilon_ + pow(a_,3.0))/max(sqrtS2*k_/epsilon_, dimensionedScalar("1.87",dimensionSet(0,0,0,0,0), 1.87))
    );

    dimensionedScalar C4Star = 2.0/Cmu_*(1.0 - C4_);
    dimensionedScalar C5Star = 2.0/Cmu_*(1.0 - C5_);
    rotSource_ = ( (C4Star*((b & S) && S) - ( C5Star *((b & W) && S)))
      /(k_*S2/epsilon_ + dimensionedScalar("small", dimensionSet(0,0,-1,0,0), VSMALL)));

    fh2_.storePrevIter();
    fh2_ =
        fmu/Ts/Cmu_*(2.0/3.0 - C3_/2.0)
        + rotSource_;

    scalar relaxFactor;
    if (this->mesh().solution().relaxField("v2kSources"))
    {
        relaxFactor = this->mesh().solution().fieldRelaxationFactor("v2kSources");
    }
    else
    {
        relaxFactor = 1.;
    }

    fh2_.relax(relaxFactor);

    if (boundV2kSource_)
    {
        scalar maxFh2 = 0.;
        forAll(fh2_, ci)
        {
            fh2_.primitiveFieldRef()[ci] = max(maxFh2, fh2_.primitiveFieldRef()[ci]);
        }
        forAll(fh2_.boundaryField(), pI)
        {
            fh2_.boundaryFieldRef()[pI].fvPatchScalarField::operator=
            (
                max(fh2_.boundaryField()[pI], maxFh2)
            );
        }
    }

    volScalarField sourceExp
    (
        alpha*rho*(pow(a_,3.0)*fh2_)
    );

    const volScalarField dfw
    (
        "dfw",
        -(Ceps2_ - 1.0 + 5.0 - 1.0/Cmu_)/Ts
    );
    volScalarField sourceImpl
    (
        alpha*rho*
        (
            ((1.0 - pow(a_,3.0))*dfw)
          - ((2.0 - Ceps1_)*nutOverk_*S2)
          - (pow(a_,3.0)*1./Ts*(C1_ + Ceps2_ -  2.0 + C2_*Sk/epsilon_))
          + (pow(a_,3.0)*C3star_/sqrt(2.0)*sqrtS2)
        )
    );

//    if (this->mesh_.time().outputTime())
//    {
//        Wshur.write();
//        S.write();
//        fmu.write();
//        rotSource_.write();
//        dfw.write();
//        G.write();
//    }

    tmp<fvScalarMatrix> v2kEqn
    (
        fvm::ddt(alpha, rho, v2k_)
      + fvm::div(alphaRhoPhi, v2k_)
      - fvm::laplacian(alpha*rho*DkEff(), v2k_)
        ==
        sourceExp
      + fvm::SuSp(sourceImpl, v2k_)
      + fvOptions(alpha, rho, v2k_)
    );

    v2kEqn.ref().relax();

    fvOptions.constrain(v2kEqn.ref());

    solve(v2kEqn);
    fvOptions.correct(v2k_);
    bound(v2k_, 0., v2kMax_);

    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
ellipticBlendingLagKE<BasicTurbulenceModel>::ellipticBlendingLagKE
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
    ellipticBlendingKEBase(),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.22
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.7
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.9
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.8
        )
    ),
    C3star_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3star",
            this->coeffDict_,
            0.65
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            this->coeffDict_,
            0.625
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.2
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.164
        )
    ),
    CT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT",
            this->coeffDict_,
            4.0
        )
    ),
    CtLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtLim",
            this->coeffDict_,
            1.0
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            75.0
        )
    ),
    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            this->coeffDict_,
            1.44
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.9   // for Lag
        )
    ),
    Ceps3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps3",
            this->coeffDict_,
            4.6
        )
    ),
    Ceps4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps4",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.2   // 1.15 in Sylvain & Flavien's paper
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
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
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    v2k_
    (
        IOobject
        (
            IOobject::groupName("v2k", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    a_
    (
        IOobject
        (
            IOobject::groupName("a", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    rotSource_
    (
        IOobject
        (
            IOobject::groupName("rotSource", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("scalar", dimless/dimTime,  0.)
    ),
    nutOverk_
    (
        IOobject
        (
            IOobject::groupName("nutOverk", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("scalar", dimTime,  0.)
    ),
    fh2_
    (
        IOobject
        (
            IOobject::groupName("fh2", U.group()),
            this->runTime_.timeName(),
            this->db_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("scalar", dimless/dimTime,  0.)
    ),
    v2kMin_(dimensionedScalar("v2kMin", v2k_.dimensions(), 1e-11)),
    v2kMax_(dimensionedScalar("v2kMax", v2k_.dimensions(), 2.0/3.0)),
    aMin_(dimensionedScalar("aMin", a_.dimensions(), 0.0)),
    fixedPointIterations_
    (
        dimensioned<label>::lookupOrAddToDict
        (
            "fixedPointIter",
            this->coeffDict_,
            1
        )
    ),
    boundV2kSource_
    (
        Switch::lookupOrAddToDict
        (
            "boundV2kSource",
            this->coeffDict_,
            false
        )
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(v2k_, v2kMin_, v2kMax_);
    bound(a_, aMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool ellipticBlendingLagKE<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C3star_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Ceps3_.readIfPresent(this->coeffDict());
        Ceps4_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        fixedPointIterations_.readIfPresent(this->coeffDict());
        boundV2kSource_.readIfPresent("boundV2kSource",this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void ellipticBlendingLagKE<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    for (int i=0; i<fixedPointIterations_.value(); ++i)
    {
        epsilon_.storePrevIter();
        k_.storePrevIter();
        a_.storePrevIter();
        v2k_.storePrevIter();

        tmp<volTensorField> tgradU = fvc::grad(this->U());

        solveTurb(tgradU());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
