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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "EdwardsSA.H"
#include "cfdTools/general/bound/bound.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> EdwardsSA<BasicTurbulenceModel>::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> EdwardsSA<BasicTurbulenceModel>::fw
(
    const volScalarField& Stilda
) const
{
    const volScalarField& y(wallDist::New(this->mesh_).y());
    volScalarField r
    (
        tanh
        (
            this->nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), ROOTVSMALL)
               )
              *sqr(this->kappa_*y)
            )
        )/tanh(1.0)
    );
    r.boundaryFieldRef().forceAssign(0.0);

    const volScalarField g(r + this->Cw2_*(pow6(r) - r));

    return g*pow((1.0 + pow6(this->Cw3_))/(pow6(g) + pow6(this->Cw3_)), 1.0/6.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
EdwardsSA<BasicTurbulenceModel>::EdwardsSA
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
    SpalartAllmaras<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void EdwardsSA<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    volScalarField& nuTilda = this->nuTilda_;

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    //prevent zero chi for Stilda calc
    bound(chi, dimensionedScalar("chimin", chi.dimensions(), SMALL), false);

    autoPtr<volScalarField> Stilda;
    {
        tmp<volTensorField> gradU(this->gradUnwc());

        Stilda.reset
        (
            new volScalarField
            (
                sqrt
                (
                    ((2*symm(gradU())) && gradU())
                    - 2.0/3.0 * sqr(tr(gradU()))
                )
                *
                (1/chi + fv1)
            )
        );
    }

    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda)
      + fvm::div(alphaRhoPhi, nuTilda)
      - fvm::laplacian(alpha*rho*this->DnuTildaEff(), nuTilda)
      - this->Cb2_/this->sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda))
     ==
        this->Cb1_*alpha*rho*Stilda()*nuTilda
      - fvm::Sp(this->Cw1_*alpha*rho*fw(Stilda())*nuTilda/sqr(y), nuTilda)
    );

    nuTildaEqn.ref().relax();
    solve(nuTildaEqn);
    bound(nuTilda, dimensionedScalar("0", nuTilda.dimensions(), 0.0));
    nuTilda.correctBoundaryConditions();

    this->correctNut(fv1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
