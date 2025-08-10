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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "functionTIso.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(functionTIso, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        functionTIso,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionTIso::functionTIso
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
    functionTIso::read();
}


Foam::autoPtr<Foam::functionTIso>
Foam::functionTIso::clone() const
{
    return autoPtr<functionTIso>
    (
        new functionTIso(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::functionTIso::castScalarModel
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


Foam::baseModels<Foam::vector>* Foam::functionTIso::castVectorModel
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


Foam::vector Foam::functionTIso::vKappaCell(const label celli) const
{
    const scalar kap(kappa(T_->operator[](celli)));
    return vector(kap, kap, kap);
}


Foam::tmp<Foam::vectorField> Foam::functionTIso::vKappaPatch
(
    const label patchi
) const
{
    tmp<vectorField> tvKappa
    (
        new vectorField(mesh_.boundaryMesh()[patchi].size(), Zero)
    );
    vectorField& vKappa = tvKappa.ref();
    tmp<scalarField> T(T_->boundaryField()[patchi]);
    forAll(vKappa, facei)
    {
        const scalar kap(kappa(T.ref()[facei]));
        vKappa[facei] = vector(kap, kap, kap);
    }
    return tvKappa;
}


Foam::tmp<Foam::vectorField> Foam::functionTIso::vKappaInternal() const
{
    tmp<vectorField> tvKappa
    (
        new vectorField(mesh_.nCells(), Zero)
    );
    vectorField& vKappa = tvKappa.ref();
    tmp<scalarField> T(T_->primitiveField());
    forAll(vKappa, celli)
    {
        const scalar kap(kappa(T.ref()[celli]));
        vKappa[celli] = vector(kap, kap, kap);
    }
    return tvKappa;
}


Foam::tmp<Foam::volVectorField>
Foam::functionTIso::vKappaGeometric() const
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


bool Foam::functionTIso::read()
{
    T_ = constructOrReturnRefFieldPtr<scalar>("T");
    kappa_.reset(Function1<scalar>::New("kappa", *dict_));

    return true;
}


// ************************************************************************* //
