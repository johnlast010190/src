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
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "RAS/buoyantKEpsilon/buoyantKEpsilon.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
buoyantKEpsilon<BasicTurbulenceModel>::buoyantKEpsilon
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
    kEpsilon<BasicTurbulenceModel>
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
            1.0
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
bool buoyantKEpsilon<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
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
buoyantKEpsilon<BasicTurbulenceModel>::Gcoef
(
    const bool deferred
) const
{
    const uniformDimensionedVectorField& g =
        this->db_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (deferred)
    {
        return
            (Cg_*this->Cmu_)*this->alpha_*this->k_*(g & fvc::grad(this->rho_))
           /(this->epsilon_.prevIter() + this->epsilonMin_);
    }
    else
    {
        return
            (Cg_*this->Cmu_)*this->alpha_*this->k_*(g & fvc::grad(this->rho_))
           /(this->epsilon_ + this->epsilonMin_);
    }
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKEpsilon<BasicTurbulenceModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->db_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > SMALL)
    {
        return -fvm::SuSp(Gcoef(this->deferredS2_), this->k_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::kSource();
    }
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->db_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > SMALL)
    {
        vector gHat(g.value()/mag(g.value()));

        volScalarField v(gHat & this->U_);
        volScalarField u
        (
            mag(this->U_ - gHat*v)
          + dimensionedScalar("SMALL", dimVelocity, SMALL)
        );

        return -fvm::SuSp(this->C1_*tanh(mag(v)/u)*Gcoef(), this->epsilon_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::epsilonSource();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
