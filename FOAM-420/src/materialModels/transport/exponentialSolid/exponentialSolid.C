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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "exponentialSolid.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponentialSolid, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        exponentialSolid,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::exponentialSolid::exponentialSolid
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    exponentialSolid::read();
}


Foam::autoPtr<Foam::exponentialSolid>
Foam::exponentialSolid::clone() const
{
    return autoPtr<exponentialSolid>
    (
        new exponentialSolid(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::exponentialSolid::castScalarModel
(
    const word& modelName
)
{
    if (modelName == kappaModel::typeName)
    {
        return dynamic_cast<kappaModel*>(this);
    }
    return nullptr;
}


Foam::baseModels<Foam::vector>* Foam::exponentialSolid::castVectorModel
(
    const word& modelName
)
{
    if (modelName == vKappaModel::typeName)
    {
        return dynamic_cast<vKappaModel*>(this);
    }
    return nullptr;
}


Foam::scalar Foam::exponentialSolid::kappaCell(const label celli) const
{
    return (kappa0_*pow(T_->operator[](celli)/Tref_, n0_));
}


Foam::tmp<Foam::scalarField> Foam::exponentialSolid::kappaPatch
(
    const label patchi
) const
{
    return (kappa0_*pow(T_->boundaryField()[patchi]/Tref_, n0_));
}


Foam::tmp<Foam::scalarField> Foam::exponentialSolid::kappaInternal() const
{
    return (kappa0_*pow(T_->primitiveField()/Tref_, n0_));
}


Foam::tmp<Foam::volScalarField>
Foam::exponentialSolid::kappaGeometric() const
{
    tmp<volScalarField> tkappa
    (
        new volScalarField
        (
            IOobject
            (
                "kappa",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("kappa", kappaModel::modelDims, Zero)
        )
    );
    volScalarField& kappa = tkappa.ref();
    kappa.primitiveFieldRef() = kappaInternal();

    forAll(kappa.boundaryField(), patchi)
    {
        kappa.boundaryFieldRef()[patchi].forceAssign(kappaPatch(patchi));
    }

    return tkappa;
}


Foam::vector Foam::exponentialSolid::vKappaCell(const label celli) const
{
    const scalar kappa(kappa0_*pow(T_->operator[](celli)/Tref_, n0_));
    return vector(kappa, kappa, kappa);
}


Foam::tmp<Foam::vectorField> Foam::exponentialSolid::vKappaPatch
(
    const label patchi
) const
{
    const scalarField kappa(kappa0_*pow(T_->boundaryField()[patchi]/Tref_, n0_));
    tmp<vectorField> tvKappa
    (
        new vectorField(kappa.size(), Zero)
    );
    vectorField& vKappa = tvKappa.ref();
    forAll(vKappa, facei)
    {
        vKappa[facei] = vector(kappa[facei], kappa[facei], kappa[facei]);
    }
    return tvKappa;
}


Foam::tmp<Foam::vectorField> Foam::exponentialSolid::vKappaInternal() const
{
    const scalarField kappa(kappa0_*pow(T_->primitiveField()/Tref_, n0_));
    tmp<vectorField> tvKappa
    (
        new vectorField(kappa.size(), Zero)
    );
    vectorField& vKappa = tvKappa.ref();
    forAll(vKappa, celli)
    {
        vKappa[celli] = vector(kappa[celli], kappa[celli], kappa[celli]);
    }
    return tvKappa;
}


Foam::tmp<Foam::volVectorField>
Foam::exponentialSolid::vKappaGeometric() const
{
    tmp<volVectorField> tvKappa
    (
        new volVectorField
        (
            IOobject
            (
                "vKappa",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedVector("vKappa", vKappaModel::modelDims, Zero)
        )
    );
    volVectorField& vKappa = tvKappa.ref();
    vKappa.primitiveFieldRef() = vKappaInternal();

    forAll(vKappa.boundaryField(), patchi)
    {
        vKappa.boundaryFieldRef()[patchi].forceAssign(vKappaPatch(patchi));
    }

    return tvKappa;
}

bool Foam::exponentialSolid::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    kappa0_ = dict_->lookup<scalar>("kappa0");
    n0_ = dict_->lookup<scalar>("n0");
    Tref_ = dict_->lookup<scalar>("Tref");

    return true;
}


// ************************************************************************* //
