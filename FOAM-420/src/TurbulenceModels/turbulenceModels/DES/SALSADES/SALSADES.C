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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "SALSADES.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SALSADES<BasicTurbulenceModel>::Cb1
(
    const volScalarField& SStar
) const
{
    const volScalarField& y(wallDist::New(this->mesh_).y());

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cb1",
                this->mesh_.time().timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            0.1335*sqrt
            (
                min
                (
                    scalar(1.25),
                    max
                    (
                        scalar(0.75),
                        max
                        (
                            pow
                            (
                                1.01*this->nuTilda_
                               /(SStar*sqr(this->kappa_*y)),
                                0.65
                            ),
                            pow(max(scalar(0.0), 1-tanh(this->chi()/68)), 0.65)
                        )
                    )
                )
            )
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SALSADES<BasicTurbulenceModel>::Omega(const volTensorField& gradU) const
{
    // NOTE: This is actually referred to as S^* in the literature, but it
    // takes the place of Omega in the standard SA model
    return sqrt(2.0)*mag(dev(symm(gradU)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SALSADES<BasicTurbulenceModel>::STilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volScalarField& SStar,
    const volScalarField& dTilda
) const
{
    return SStar * (1.0/chi + this->fv1(chi));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SALSADES<BasicTurbulenceModel>::r
(
    const volScalarField& visc,
    const volScalarField& Stilda,
    const volScalarField& dTilda
) const
{
    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<volScalarField> r
    (
        0.7*visc
       /(
           max
           (
               Stilda,
               dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
           )
           *sqr(this->kappa_*y)
         + dimensionedScalar
           (
               "ROOTVSMALL",
               dimensionSet(0, 2 , -1, 0, 0),
               ROOTVSMALL
           )
        )
    );

    if (rho0_.valid())
    {
        r = 1.6*tanh(sqrt(rho0_/this->rho())*r);
    }
    else
    {
        r = 1.6*tanh(r);
    }

    r->boundaryFieldRef().forceAssign(0.0);

    return r;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SALSADES<BasicTurbulenceModel>::SALSADES
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
    SpalartAllmarasDES<BasicTurbulenceModel>
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
    rho0_()
{
    this->sigmaNu_ =
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNu",
            this->coeffDict_,
            1.0
        )
    );

    if (this->coeffDict_.found("rho0"))
    {
        rho0_.reset
        (
            new dimensioned<scalar>("rho0", dimDensity, this->coeffDict_)
        );
    }

    this->lowReCorrection_ =
    (
        Switch::lookupOrAddToDict
        (
            "lowReCorrection",
            this->coeffDict_,
            false
        )
    );

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SALSADES<BasicTurbulenceModel>::read()
{
    if (SpalartAllmarasDES<BasicTurbulenceModel>::read())
    {
        if (this->coeffDict_.found("rho0"))
        {
            rho0_.reset
            (
                new dimensioned<scalar>
                (
                    "rho0", dimDensity, this->coeffDict_
                )
            );
        }

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
