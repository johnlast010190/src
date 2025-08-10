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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "kEpsilon.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/bound/bound.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kEpsilon<BasicTurbulenceModel>::correctNut()
{
    if (limitNut_)
    {
    // update nut according to kOmegaSSTBase implementation
    // which stabilizes nut for low-Re

        // k-omega-SST coeffs
        scalar a1(0.31);
        scalar b1(1.0);
        //scalar betaStar(0.09); //equal to Cmu_!!!

        // compute omega from epsilon
        volScalarField omega("omega", epsilon_/(Cmu_*k_));

        // compute S3 and F23
        tmp<volScalarField> S2 = 2*magSqr(symm(fvc::grad(this->U_)));
        const volScalarField& y(wallDist::New(this->mesh_).y());
        tmp<volScalarField> arg2 =
            max
            (
                (scalar(2)/Cmu_)*sqrt(k_)
               /max(omega*y, this->omegaMin_*this->ySmall_),
                scalar(500)*(this->mu()/this->rho_)
               /max(sqr(y)*omega, sqr(this->ySmall_)*this->omegaMin_)
            );
        tmp<volScalarField> F23 = tanh(sqr(arg2));
        tmp<volScalarField> epsLim = Cmu_*b1/a1*k_*F23()*sqrt(S2);
        dimensionedScalar epsLimMax("e", epsilon_.dimensions(), epsLimMax_);

        // Correct the turbulence viscosity
        this->nut_ = Cmu_*sqr(k_)/max(epsilon_,min(epsLim,epsLimMax));
    }
    else
    {
    // standard nut update for k-epsilon
        this->nut_ = Cmu_*sqr(k_)/epsilon_;
    }

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_, this->nut_.db()).correct(this->nut_);

    BasicTurbulenceModel::correctNut();

    //- bound nut. Can be triggered in bad quality meshes in the beginning of
    // the simulation
    scalar maxNu = gMax(this->nu()().primitiveField());
    dimensionedScalar nutMax
    (
        "nutMax",
        this->nut_.dimensions(),
        this->maxTurbViscRatio_*maxNu
    );
    dimensionedScalar nutMin
    (
        "nutMin",
        this->nut_.dimensions(),
        0.0
    );
    bound(this->nut_, nutMin, nutMax, false);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilon<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilon<BasicTurbulenceModel>::kEpsilon
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

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
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
            1.3
        )
    ),
    maxTurbViscRatio_
    (
        this->coeffDict_.lookupOrAddDefault
        (
            "maxTurbViscosityRatio",
            1e+6
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
    deferredS1_
    (
        Switch::lookupOrAddToDict
        (
            "deferredS1",
            this->coeffDict_,
            false
        )
    ),
    deferredS2_
    (
        Switch::lookupOrAddToDict
        (
            "deferredS2",
            this->coeffDict_,
            false
        )
    ),
    limitNut_
    (
        Switch::lookupOrAddToDict
        (
            "limitNut",
            this->coeffDict_,
            false
        )
    ),
    epsLimMax_
    (
        this->coeffDict_.lookupOrAddDefault
        (
            "epsLimMax",
            1e-7
        )
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    if (this->coeffDict_.found("deferred"))
    {
        bool deferred(this->coeffDict_.lookup("deferred"));
        if (deferred)
        {
            deferredS1_ = true;
            deferredS2_ = true;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilon<BasicTurbulenceModel>::L() const
{
    scalar Cmu_75(::pow(Cmu_.value(), 0.75));

    return
    (
        k_*sqrt(k_)/epsilon_*Cmu_75
    );
}

template<class BasicTurbulenceModel>
bool kEpsilon<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        this->coeffDict().readIfPresent
        (
            "maxTurbViscosityRatio",
            maxTurbViscRatio_
        );
        deferredS1_.readIfPresent("deferredS1", this->coeffDict());
        deferredS2_.readIfPresent("deferredS2", this->coeffDict());
        limitNut_.readIfPresent("limitNut", this->coeffDict());
        this->coeffDict().readIfPresent
        (
            "epsLimMax",
            epsLimMax_
        );

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kEpsilon<BasicTurbulenceModel>::correct()
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

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // store epsilon prevIter for time consistent discretisation
    if (deferredS1_ || deferredS2_)
    {
        if (debug)
        {
            Info<< "using time consistent turbulence source discretisation" << endl;
        }
        epsilon_.storePrevIter();
    }

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(deferredS1_), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
kEpsilon<BasicTurbulenceModel>::epsilonByk
(
    const bool deferred
) const
{
    if (deferred)
    {
        return epsilon_.prevIter()/k_();
    }
    else
    {
        return epsilon_()/k_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
