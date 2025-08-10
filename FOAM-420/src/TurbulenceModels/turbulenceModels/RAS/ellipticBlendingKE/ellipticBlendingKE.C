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

#include "ellipticBlendingKE.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/bound/bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> ellipticBlendingKE<BasicTurbulenceModel>::Ts() const
{
    //  ( Teddy + Tkol )
    return sqrt((k_*k_/sqr(epsilon_ + dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL))) +
                CT_*CT_*(max(this->nu(), dimensionedScalar("zero", this->nu()().dimensions(), 0.0))
               /(epsilon_ + dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL))));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> ellipticBlendingKE<BasicTurbulenceModel>::TsLim() const
{
    //  ( Treal )
    return CtLim_/(sqrt(6*magSqr(dev(symm(fvc::grad(this->U_)))))*Cmu_*v2k_ +
               dimensionedScalar("TinvMin", dimensionSet(0,0,-1,0,0,0,0), SMALL));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> ellipticBlendingKE<BasicTurbulenceModel>::Ls() const
{
//    return CL_*sqrt((pow(k_,3.0)/pow(epsilon_ + dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL), 2.0)) +
//               Ceta_*Ceta_*sqrt(pow(this->nu(),3.0)/(epsilon_ +
//               dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL))));
    return CL_*sqrt((pow(k_,3.0)/pow(epsilon_ + dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL), 2.0)) +
               Ceta_*Ceta_*pow(pow(this->nu(),3.0)/(epsilon_ +
               dimensionedScalar("epsilon0", epsilon_.dimensions(), SMALL)), 0.5));
}


template<class BasicTurbulenceModel>
void ellipticBlendingKE<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2
)
{
    this->nut_ = min(this->Cmu_*k_*v2k_*Ts(), this->Cmu_*k_*v2k_*TsLim());
//    this->nut_ = this->Cmu_*v2k_*k_*Ts();
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_, this->nut_.db()).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void ellipticBlendingKE<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(dev(symm(fvc::grad(this->U_)))));
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
ellipticBlendingKE<BasicTurbulenceModel>::ellipticBlendingKE
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
//            0.6
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
            1.83
        )
    ),
    Ceps3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps3",
            this->coeffDict_,
            0.4
//            1.0
        )
    ),
    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            2.3
//            4.6
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
            1.5
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
    v2kMin_(dimensionedScalar("v2kMin", v2k_.dimensions(), SMALL)),
    v2kMax_(dimensionedScalar("v2kMax", v2k_.dimensions(), 2.0/3.0)),
    aMin_(dimensionedScalar("aMin", a_.dimensions(), 0.0))
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(v2k_, v2kMin_);
    bound(a_, aMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool ellipticBlendingKE<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Ceps3_.readIfPresent(this->coeffDict());
        Ck_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void ellipticBlendingKE<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_, U.db()));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(2*magSqr(dev(symm(gradU))));

    const volScalarField G(this->GName(), nut*S2);
    const volScalarField Ts(this->Ts());
    const volScalarField L2(type() + ":L2", sqr(Ls()));

    const volScalarField dfw
    (
       -epsilon_/2.0/k_
//            -k_/2.0/epsilon_
    );

    const volScalarField fh
    (
       -1.0/Ts*(C1_ - 1.0 + C2_*G/epsilon_)*(v2k_ - 2.0/3.0) //Not matching the paper
//       -1.0/Ts*(C1_ - 1.0 + C2_*G/k_)*(v2k_ - 2.0/3.0)
    );
    const volScalarField Ceps2star
    (
        "Ceps2star",
        Ceps2_ + pow(a_,3)*(Ceps3_ - Ceps2_)*tanh(pow(mag(fvc::div(nut/sigmaK_*fvc::grad(k_))/epsilon_),1.5))
    );
    const volSymmTensorField S(symm(gradU));
//    const volVectorField& n_(wallDist::New(this->mesh_).n());
    const volVectorField grada(fvc::grad(a_));
    const volVectorField n(grada/(mag(grada) + dimensionedScalar("nsmall", grada.dimensions(), VSMALL)));

//    const volScalarField E
//    (
//        2.0*this->nu()*nut*sqr(fvc::div(mag(2*S & n) * n )
//        )
//    );
//
    const volScalarField E
    (
       2.0*this->nu()*nut*fvc::magSqrGradGrad(U)
    );

    const volScalarField Se
    (
        "Se",
        Ceps1_*alpha*rho*G/Ts
    );
    const volScalarField eCoeff
    (
        "eCoeff",
        Ceps2star*alpha*rho/Ts
    );
    // Update epsilon (and possibly G) at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
//        Ceps1_*alpha*rho*G*epsilon_/k //starCCM
//        - fvm::Sp(Ceps2star*alpha*rho*epsilon_/k_, epsilon_) //starCCM
        Se
//        Ceps1_*alpha*rho*G/Ts   // converges better with this option
//      - fvm::Sp(Ceps2star*alpha*rho/Ts, epsilon_)
      - fvm::Sp(eCoeff, epsilon_)
//      + Ck_*pow(1.0 - a_,3.0)*E          // E term in epsilon equation
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    scalar unityCoef = 1;

    const volScalarField Sk
    (
        "Sk",
        unityCoef*G
    );
    volScalarField ECoeff
    (
        "ECoeff",
        Ck_*pow(1.0 - a_,3.0)*E/epsilon_
    );
    volScalarField EpsCoeff
    (
        "EpsCoeff",
        epsilon_/k_
    );

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*Sk
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_) // Unknown term??
      - fvm::Sp(alpha*rho*EpsCoeff, k_)
      - fvm::Sp(ECoeff, k_)  // try without
      + fvOptions(alpha, rho, k_)
    );


    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    // Elliptic Blending function equation (a is used because alpha is taken by void fraction)
    tmp<fvScalarMatrix> aEqn
    (
      - fvm::laplacian(a_)
      - 1.0/L2
     ==
      - fvm::Sp(1.0/L2, a_)
    );

    aEqn.ref().relax();
    fvOptions.constrain(aEqn.ref());
    solve(aEqn);
    fvOptions.correct(a_);
    bound(a_, aMin_);

    // Turbulence stress normal to streamlines equation
    tmp<fvScalarMatrix> v2kEqn
    (
        fvm::ddt(alpha, rho, v2k_)
      + fvm::div(alphaRhoPhi, v2k_)
      - fvm::laplacian(alpha*rho*DkEff(), v2k_)
      ==
        alpha*rho*(pow(a_,3.0)*fh)
      + fvm::Sp(alpha*rho*((1.0 - pow(a_,3.0))*dfw), v2k_)
      - fvm::Sp(alpha*rho*G/k_, v2k_)
//      + 2.0/k_*nut/sigmaK_*(fvc::grad(a_) & fvc::grad(k_)) //not matching the paper
      + 2.0/k_*nut/sigmaK_*(fvc::grad(k_) & fvc::grad(v2k_))
      + fvOptions(alpha, rho, v2k_)
    );

    v2kEqn.ref().relax();
    fvOptions.constrain(v2kEqn.ref());
    solve(v2kEqn);
    fvOptions.correct(v2k_);
    bound(v2k_, v2kMin_);
    if (max(v2k_).value() > v2kMax_.value())
    {
        v2k_.primitiveFieldRef() = min
        (
            min
            (
                v2k_.primitiveField(),
                fvc::average(min(v2k_, v2kMax_))().primitiveField()
               * neg(-v2k_.primitiveField())
            ),
            v2kMax_.value()
        );
        Info<< "new max: " << gMax(v2k_.primitiveField()) << endl;
        v2k_.boundaryFieldRef() = min(v2k_.boundaryField(), v2kMax_.value());
    }
    v2k_ = min(v2k_ ,v2kMax_);
    correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
