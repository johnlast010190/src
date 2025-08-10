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
    (c) 2010-2016 Esi Ltd.
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "linearViscousStress.H"
#include "finiteVolume/fvc/fvc.H"
#include "finiteVolume/fvm/fvm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::linearViscousStress<BasicTurbulenceModel>::linearViscousStress
(
    const word& modelName,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        modelName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::linearViscousStress<BasicTurbulenceModel>::read()
{
    return BasicTurbulenceModel::read();
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::linearViscousStress<BasicTurbulenceModel>::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->U_.group()),
                this->runTime_.timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (-(this->alpha_*this->rho_*this->nuEff()))
           *dev(twoSymm(fvc::grad(this->U_)))
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    U.updateBoundaryCoeffs();
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U,
    const word& scheme
) const
{
    U.updateBoundaryCoeffs();
    return
    (
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
      - fvc::div
      (
          (this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))),
          scheme
      )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    U.updateBoundaryCoeffs();
    return
    (
      - fvc::div((this->alpha_*rho*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U,
    const word& scheme
) const
{
    U.updateBoundaryCoeffs();
    return
    (
      - fvc::div
        (
            (this->alpha_*rho*this->nuEff())*dev2(T(fvc::grad(U))),
            scheme
        )
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::linearViscousStress<BasicTurbulenceModel>::bDivDevRhoReff
(
    volVectorField& U
) const
{
    U.updateBoundaryCoeffs();
    return
    (
      - fvm::bLaplacian(this->alpha_*this->rho_*this->nuEff(), U)
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector>>
Foam::linearViscousStress<BasicTurbulenceModel>::bDivDevRhoReff
(
    volVectorField& U,
    const word& scheme
) const
{
    U.updateBoundaryCoeffs();
    return
    (
      - fvm::bLaplacian(this->alpha_*this->rho_*this->nuEff(), U)
      - fvc::div
      (
          (this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))),
          scheme
      )
    );
}


template<class BasicTurbulenceModel>
void Foam::linearViscousStress<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
