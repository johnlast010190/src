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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "kOmegaSSTDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDDES<BasicTurbulenceModel>::rd
(
    const volScalarField& magGradU
) const
{
    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<volScalarField> tr
    (
        min
        (
            this->nuEff()
           /(
                max
                (
                    magGradU,
                    dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
                )
                *sqr(this->kappa_*y)
            ),
            scalar(10)
        )
    );
    tr.ref().boundaryFieldRef().forceAssign(0.0);

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow(Cd1_*rd(magGradU), Cd2_));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    const volScalarField lRAS(sqrt(k)/(this->betaStar_*omega));
    const volScalarField lLES(CDES*this->delta());

    return max
    (
        lRAS
      - fd(magGradU)
       *max
        (
            lRAS - lLES,
            dimensionedScalar("zero", dimLength, 0)
        ),
        dimensionedScalar("small", dimLength, SMALL)
    );
}

template<class BasicTurbulenceModel>
void kOmegaSSTDDES<BasicTurbulenceModel>::correct()
{
        // Correct the turbulence viscosity
    kOmegaSSTDES<BasicTurbulenceModel>::correct();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTDDES<BasicTurbulenceModel>::kOmegaSSTDDES
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

    Cd1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd1",
            this->coeffDict_,
            20
        )
    ),
    Cd2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd2",
            this->coeffDict_,
            3
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
bool kOmegaSSTDDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTDES<BasicTurbulenceModel>::read())
    {
        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());

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
