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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTSASDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTSASDES<BasicTurbulenceModel>::Lvk2
(
    const volScalarField& S2
) const
{
    return max
    (
        this->kappa_*sqrt(S2)
       /(
            mag(fvc::laplacian(this->U()))
          + dimensionedScalar
            (
                "ROOTVSMALL",
                dimensionSet(0, -1 , -1, 0, 0, 0, 0),
                ROOTVSMALL
            )
        ),
        this->Cs_*this->delta()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTSASDES<BasicTurbulenceModel>::kOmegaSSTSASDES
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
    kOmegaSSTDES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.262
        )
    ),
    alphaPhi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaPhi",
            this->coeffDict_,
            0.666667
        )
    ),
    zetaTilda2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "zetaTilda2",
            this->coeffDict_,
            1.755
        )
    ),
    FSAS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "FSAS",
            this->coeffDict_,
            1.25
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTSASDES<BasicTurbulenceModel>::Qsas
(
    const volScalarField& S2,
    const volScalarField& gamma,
    const volScalarField& beta
) const
{
    volScalarField L(sqrt(this->k_)/(pow025(this->betaStar_)*this->omega_));
    tmp<volScalarField> grad_omega_k = max
    (
        magSqr(fvc::grad(this->omega_))/sqr(this->omega_),
        magSqr(fvc::grad(this->k_))/sqr(this->k_)
    );

    return fvm::Su
        (
            this->alpha_*this->rho_*FSAS_
           *max
            (
                dimensionedScalar("zero",dimensionSet(0, 0, -2, 0, 0), 0.0),
                zetaTilda2_*this->kappa_*S2*sqr(L/Lvk2(S2))
              - 2.0/alphaPhi_*this->k_*grad_omega_k
            ),
            this->omega_
        );
}


template<class BasicTurbulenceModel>
bool kOmegaSSTSASDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTDES<BasicTurbulenceModel>::read())
    {
        Cs_.readIfPresent(this->coeffDict());
        alphaPhi_.readIfPresent(this->coeffDict());
        zetaTilda2_.readIfPresent(this->coeffDict());
        FSAS_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
