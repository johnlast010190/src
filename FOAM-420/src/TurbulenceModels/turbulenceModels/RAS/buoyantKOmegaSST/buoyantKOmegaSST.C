/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM® : Professional Open-source CFD
|   o   O   o    |  Version : 4.2.0
|    o     o     |  Copyright © 2019 ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM® <http://www.openfoam.org/>.

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

Portions Copyright © 2014-2015 OpenFOAM Foundation.

Author
    2018. Lisandro Maders (Esi Ltd.). All rights reserved.

\*---------------------------------------------------------------------------*/

#include "buoyantKOmegaSST.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
buoyantKOmegaSST<BasicTurbulenceModel>::buoyantKOmegaSST
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
    kOmegaSST<BasicTurbulenceModel>
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
    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1.1765
        )
    ),
    C1o_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1o",
            this->coeffDict_,
            1.0
        )
    ),
    C2o_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2o",
            this->coeffDict_,
            1.0
        )
    ),
    Gb_
    (
        IOobject
        (
            "Gb",
            this->runTime_.timeName(),
            this->db_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Gcoef()
     ),
    diss_
    (
        Switch::lookupOrAddToDict
        (
            "dissipation",
            this->coeffDict_,
            true
        )
    )
{

    if (phi.dimensions() == dimMass/dimTime)
    {
        compressible_ = true;
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        compressible_ = false;
    }

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool buoyantKOmegaSST<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        Cg_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
buoyantKOmegaSST<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
    this->db_.objectRegistry::template
    lookupObject<uniformDimensionedVectorField>("g");

    return (-1*this->nut_*Cg_*(g & fvc::grad(this->rho()))) ;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKOmegaSST<BasicTurbulenceModel>::kSource() const
{

    const uniformDimensionedVectorField& g =
        this->db_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > SMALL)
    {
        dimensionedScalar kSmall
        (
            "kSmall",
            this->k_.dimensions(),
            SMALL
        );

        if (compressible_)
        {
            return fvm::Sp(Gcoef()/max(this->k_, kSmall), this->k_);
        }
        else
        {
            return fvm::Sp(Gcoef()/this->rho()/max(this->k_, kSmall), this->k_);
        }
    }
    else
    {
        return kOmegaSST<BasicTurbulenceModel>::kSource();
    }

}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKOmegaSST<BasicTurbulenceModel>::omegaSource() const
{
    const uniformDimensionedVectorField& g =
        this->db_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > SMALL && diss_)
    {
        dimensionedScalar cone
        (
            "cone",
            this->gamma1_.dimensions(),
            1.0
        );

        dimensionedScalar czero
        (
            "czero",
            this->Gb_.dimensions(),
            0.0
        );

        vector gHat(g.value()/mag(g.value()));

        volScalarField v(gHat & this->U_);
        volScalarField u
        (
            mag(this->U_ - gHat*v)
          + dimensionedScalar("SMALL", dimVelocity, SMALL)
        );

        //calculate blended coeff. SST model
        tmp<volScalarField> CDkOmega =
        (
                (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
        );
        tmp<volScalarField::Internal>bgamma(this->gamma(this->F1(CDkOmega())));

        tmp<volScalarField>Gcf(Gcoef());

        tmp<volScalarField> klim = max
        (
            this->k_,
            dimensionedScalar("0", this->k_.dimensions(), 0.0)
        );

        tmp<volScalarField> C3e = tanh
        (
            tanh(mag(v)/u)
        );

        return fvm::SuSp
        (
            ((cone+bgamma())*C1o_*C3e->internalField()*max(Gcf->internalField(), czero)-C2o_*Gcf->internalField())
            /klim->internalField(),
            this->omega_
        );

    }
    else
    {
        return kOmegaSST<BasicTurbulenceModel>::omegaSource();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
