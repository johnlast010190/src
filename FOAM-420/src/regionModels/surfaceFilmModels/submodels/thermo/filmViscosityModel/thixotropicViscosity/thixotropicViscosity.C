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
    (c) 2013-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/thermo/filmViscosityModel/thixotropicViscosity/thixotropicViscosity.H"
#include "kinematicSingleLayer/kinematicSingleLayer.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "finiteVolume/fvm/fvmDdt.H"
#include "finiteVolume/fvm/fvmDiv.H"
#include "finiteVolume/fvc/fvcDiv.H"
#include "finiteVolume/fvm/fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thixotropicViscosity, 0);

addToRunTimeSelectionTable
(
    filmViscosityModel,
    thixotropicViscosity,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thixotropicViscosity::thixotropicViscosity
(
    surfaceFilmModel& film,
    const dictionary& dict,
    volScalarField& mu
)
:
    filmViscosityModel(typeName, film, dict, mu),
    a_("a", dimless/dimTime, coeffDict_),
    b_("b", dimless, coeffDict_),
    d_("d", dimless, coeffDict_),
    c_("c", pow(dimTime, d_.value() - scalar(1)), coeffDict_),
    mu0_("mu0", dimPressure*dimTime, coeffDict_),
    muInf_("muInf", mu0_.dimensions(), coeffDict_),
    K_(1.0 - Foam::sqrt(muInf_/mu0_)),
    lambda_
    (
        IOobject
        (
            typeName + ":lambda",
            film.regionMesh().time().timeName(),
            film.regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh()
    )
{
    lambda_.min(1.0);
    lambda_.max(0.0);

    // Initialise viscosity to inf value because it cannot be evaluated yet
    mu_ = muInf_;
    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thixotropicViscosity::~thixotropicViscosity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thixotropicViscosity::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    const kinematicSingleLayer& film = filmType<kinematicSingleLayer>();

    // references to film fields
    const volVectorField& U = film.U();
    const volVectorField& Uw = film.Uw();
    const volScalarField& delta = film.delta();
    const volScalarField& deltaRho = film.deltaRho();
    const surfaceScalarField& phi = film.phi();
    const volScalarField& alpha = film.alpha();
    const Time& runTime = this->film().regionMesh().time();

    // Shear rate
    volScalarField gDot("gDot", alpha*mag(U - Uw)/(delta + film.deltaSmall()));

    if (debug && runTime.writeTime())
    {
        gDot.write();
    }

    dimensionedScalar deltaRho0("deltaRho0", deltaRho.dimensions(), ROOTVSMALL);
    surfaceScalarField phiU(phi/fvc::interpolate(deltaRho + deltaRho0));

    dimensionedScalar c0("c0", dimless/dimTime, ROOTVSMALL);
    volScalarField coeff("coeff", -c_*pow(gDot, d_) + c0);

    // Limit the filmMass and deltaMass to calculate the effect of the added
    // droplets with lambda = 0 to the film
    const volScalarField filmMass
    (
        "thixotropicViscosity:filmMass",
        film.netMass() + dimensionedScalar("SMALL", dimMass, ROOTVSMALL)
    );
    const volScalarField deltaMass
    (
        "thixotropicViscosity:deltaMass",
        max(dimensionedScalar("zero", dimMass, 0), film.deltaMass())
    );

    fvScalarMatrix lambdaEqn
    (
        fvm::ddt(lambda_)
      + fvm::div(phiU, lambda_)
      - fvm::Sp(fvc::div(phiU), lambda_)
      ==
        a_*pow((1.0 - lambda_), b_)
      + fvm::SuSp(coeff, lambda_)
      - fvm::Sp(deltaMass/(runTime.deltaT()*filmMass), lambda_)
    );

    lambdaEqn.relax();
    lambdaEqn.solve();

    lambda_.min(1.0);
    lambda_.max(0.0);

    mu_ = muInf_/(sqr(1.0 - K_*lambda_) + ROOTVSMALL);
    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
