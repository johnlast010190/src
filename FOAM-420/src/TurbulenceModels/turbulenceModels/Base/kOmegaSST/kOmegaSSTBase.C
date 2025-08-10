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
    (c) 2010-2017, 2020 Esi Ltd
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTBase.H"
#include "cfdTools/general/bound/bound.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTBase<BasicEddyViscosityModel>::productionTerm
(
    const tmp<volTensorField>& tgradU,
    const tmp<volTensorField>& tSkew,
    const tmp<volSymmTensorField>& tSymm
) const
{

    if (modelType_ == "standard")
    {
        return tmp<volScalarField::Internal>
               (
                   new volScalarField::Internal
                   (
                       tgradU() && dev(twoSymm(tgradU()))
                   )
               );
    }
    else if (modelType_ == "vorticity")
    {
        return tmp<volScalarField::Internal>
               (
                   new volScalarField::Internal
                   (
                       2*magSqr(tSkew())
                   )
               );
    }
    else if (modelType_ == "blended")
    {
        return tmp<volScalarField::Internal>
               (
                   new volScalarField::Internal
                   (
                       2*mag(symm(tgradU()))*mag(tSkew())
                   )
               );
    }
    else
    {
        FatalErrorInFunction << "productionType " << modelType_
            << " is not supported. " << endl
            << "Please specify either: standard | vorticity | blended"
            << exit(FatalError);

        return tmp<volScalarField::Internal>();
    }
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTBase<BasicEddyViscosityModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)
               /max(omega_*y, this->omegaMin_*this->ySmall_),
                scalar(500)*(this->mu()/this->rho_)
               /max(sqr(y)*omega_, sqr(this->ySmall_)*this->omegaMin_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTBase<BasicEddyViscosityModel>::F2() const
{
    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<volScalarField> arg2 =
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)
           /max(omega_*y, this->omegaMin_*this->ySmall_),
            scalar(500)*(this->mu()/this->rho_)
           /max(sqr(y)*omega_, sqr(this->ySmall_)*this->omegaMin_)
        );
    //NOTE: Esi removed min(...,100) function around arg2 in compressible
    // version but not incompressible version...as I don't see any mention of
    // this at http://turbmodels.larc.nasa.gov/sst.html, I have removed it here

    return tanh(sqr(arg2));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTBase<BasicEddyViscosityModel>::F3() const
{
    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTBase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    correctNut(S2, true);
}


template<class BasicEddyViscosityModel>
void kOmegaSSTBase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(true);
}


template<class BasicEddyViscosityModel>
void kOmegaSSTBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2,
    bool createFvOptions
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));

    if (!limitNearWallCell_)
    {
        forAll(this->nut_.boundaryField(), patchi)
        {
            const fvPatchScalarField& nutwpf = this->nut_.boundaryField()[patchi];
            if (isA<nutWallFunctionFvPatchScalarField>(nutwpf))
            {
                const auto& nutWF =
                    dynamic_cast<const nutWallFunctionFvPatchScalarField&>(nutwpf);
                forAll(nutWF, facei)
                {
                    label nearWallCelli = nutWF.patch().faceCells()[facei];
                    this->nut_[nearWallCelli] = k_[nearWallCelli]/omega_[nearWallCelli];
                }
            }
        }
    }
    this->nut_.correctBoundaryConditions();

    if (createFvOptions)
    {
        fv::options::New(this->mesh_, this->nut_.db()).correct(this->nut_);
    }

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


template<class BasicEddyViscosityModel>
void kOmegaSSTBase<BasicEddyViscosityModel>::correctNut(bool createFvOptions)
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), createFvOptions);
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTBase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G,
    const bool deferred
) const
{
    if (deferred)
    {
        return min(G, (c1_*betaStar_)*this->k_()*this->omega_.prevIter()());
    }
    else
    {
        return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
    }
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTBase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU,
    const bool deferred
) const
{
    if (deferred)
    {
        return betaStar_*omega_.prevIter()();
    }
    else
    {
        return betaStar_*omega_();
    }
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTBase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega_()
       *max(a1_*omega_(), b1_*F2*sqrt(S2))
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTBase<BasicEddyViscosityModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTBase<BasicEddyViscosityModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTBase<BasicEddyViscosityModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

// Based on model kOmegaSST-RC - https://turbmodels.larc.nasa.gov/sst.html#sst-rc
template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTBase<BasicEddyViscosityModel>::correctCurvature
(
    const rhoField& rho,
    const surfaceScalarField& phi,
    tmp<volTensorField> tSkew,
    tmp<volSymmTensorField> tSymm
)
{
    tmp<volScalarField> symInnerProduct(2.0*tSymm() && tSymm());
    tmp<volScalarField> asymInnerProduct
    (
        max(2.0*tSkew() && tSkew(),
        dimensionedScalar("0", dimensionSet(0, 0, -2, 0, 0), 0.0))
    );

    tmp<volScalarField> w(
        atan(dimensionedScalar(
                 "4",
                 dimensionSet(0, 0, 2, 0, 0),
                 1.0e-02) *
             asymInnerProduct()) *
            2.0 / (constant::mathematical::pi)*(asymInnerProduct() - symInnerProduct()) +
        symInnerProduct());

    tmp<volScalarField> rStar(
        sqrt(symInnerProduct() / max(
            w, dimensionedScalar("minw", w().dimensions(), SMALL))));

    tmp<volScalarField> D(sqrt(max(symInnerProduct(), 0.09*omega_*omega_)));

    tmp<volSymmTensorField> divS(fvc::DDt(rho, phi ,tSymm()));

    tmp<volScalarField> rT((tSkew().T() & tSymm()) && divS());

    divS.clear();
    tSkew.clear();
    tSymm.clear();

    tmp<volScalarField> w2(
        atan(dimensionedScalar("1",
                               dimensionSet(0, 0, 2, 0, 0), 1.0e-2) *
             asymInnerProduct()) *
            2.0 / (constant::mathematical::pi) *
            (sqrt(asymInnerProduct()) - D()) + D()
    );

    asymInnerProduct.clear();

    tmp<volScalarField> rTilda
    (
        2.0*rT()
        /(max((w2()*D()*D()*D()), dimensionedScalar("limit", dimensionSet(0, 0, -4, 0, 0), SMALL)))
    );
    symInnerProduct.clear();
    D.clear();
    w2.clear();
    rT.clear();

    tmp<volScalarField> Fr1 //curvature correction term
    (
        max
        (
            min
            (
                (scalar(1.0) + cr1_) * 2.0 * rStar() / (scalar(1) + rStar()) * (scalar(1.0)
                - cr3_ * atan(cr2_ * rTilda())) - cr1_,
                scalar(1.25)
            ),
            scalar(0.0)
        )
    );

    w.clear();
    rStar.clear();
    rTilda.clear();

    return Fr1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
kOmegaSSTBase<BasicEddyViscosityModel>::kOmegaSSTBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
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

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
    limitNearWallCell_
    (
        Switch::lookupOrAddToDict
        (
            "limitNearWallCell",
            this->coeffDict_,
            false
        )
    ),
    modelType_
    (
        this->coeffDict_.template lookupOrAddDefault<word>
        (
            "productionType",
            "standard"
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

    // curvature correction variables for kOmegaSST-RC model
    curvature_
    (
        Switch::lookupOrAddToDict
        (
            "curvature",
            this->coeffDict_,
            false
        )
    ),
    cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict(
        "cr1",
        this->coeffDict_,
        1.0)
    ),
    cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict(
        "cr2",
        this->coeffDict_,
        2.0)
    ),
    cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict(
        "cr3",
        this->coeffDict_,
        1.0)
    ),
    // --------- end of curvature correction variables -----

//    y_(wallDist::New(this->mesh_).y()),

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
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTBase<BasicEddyViscosityModel>::L() const
{
    scalar Cmu_m25(::pow(betaStar_.value(), -0.25));

    return
    (
        sqrt(k_)/omega_*Cmu_m25
    );
}


template<class BasicEddyViscosityModel>
bool kOmegaSSTBase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
        deferredS1_.readIfPresent("deferredS1", this->coeffDict());
        deferredS2_.readIfPresent("deferredS2", this->coeffDict());

        this->coeffDict().readIfPresent("productionType", modelType_);

        // --- read curvature correction related variables
        curvature_.readIfPresent("curvature", this->coeffDict());

        if (curvature_)
        {
            Info<< "Activating curvature correction of kOmegaSST model -"
                << " Smirnov, P. E., Menter, F. R. (2009)" << endl << endl;

            cr1_.readIfPresent(this->coeffDict());
            cr2_.readIfPresent(this->coeffDict());
            cr3_.readIfPresent(this->coeffDict());
        }
        // -----------------------------------------

        this->coeffDict().readIfPresent
        (
            "maxTurbViscosityRatio",
            maxTurbViscRatio_
        );

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicEddyViscosityModel>
void kOmegaSSTBase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_, U.db()));

    BasicEddyViscosityModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    tmp<volTensorField> tSkew = skew(tgradU());
    tmp<volSymmTensorField> tSymm = symm(tgradU());
    volScalarField S2(2*magSqr(tSymm()));

    tmp<volScalarField::Internal> GbyNu0t = productionTerm
    (
        tgradU,
        tSkew,
        tSymm
    );
    volScalarField::Internal& GbyNu0 = GbyNu0t.ref();

    volScalarField::Internal G(this->GName(), nut()*GbyNu0);

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // store omega prevIter for time consistent discretisation
    if (deferredS1_ || deferredS2_)
    {
        if (turbulenceModel::debug)
        {
            Info<< "using time consistent turbulence source discretisation" << endl;
        }
        omega_.storePrevIter();
    }

    tmp<volScalarField> Fr1
    (
        new volScalarField
        (
            IOobject
            (
                "Fr1",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("Fr1", dimless, scalar(1.0)),
            fixedValueFvPatchField<scalar>::typeName
        )
    );

    if (curvature_)
    {
        Fr1 = correctCurvature(rho, phi, tSkew, tSymm);
    }

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());
    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma*Fr1()*GbyNu(GbyNu0, F23(), S2()) //Fr1_ = curvature correction
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G, deferredS1_)*Fr1()             //Fr1_ = curvature correction
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU(), deferredS2_), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    tgradU.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
