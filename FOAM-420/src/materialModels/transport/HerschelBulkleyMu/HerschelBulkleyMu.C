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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "HerschelBulkleyMu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HerschelBulkleyMu, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        HerschelBulkleyMu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HerschelBulkleyMu::HerschelBulkleyMu
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
    sMod_.setSize(modelsEnumSize_);
    dep_.setSize(1);
    HerschelBulkleyMu::read();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HerschelBulkleyMu::updateTable(const word& modelName)
{
    const matScalarTable& models =
        materialTables_.sTable(phaseName_, specieName_);

    // Model names
    const word& muName = muModel::typeName;
    const word& strainRateName = strainRateModel::typeName;

    if (modelName == muName)
    {
        sMod_.set(sRateInd, models[strainRateName]);
        dep_[0].model = models[muName];
        dep_[0].dependencies.setSize(1);
        dep_[0].dependencies.set(0, models[strainRateName]);
    }
}


Foam::baseModels<Foam::scalar>* Foam::HerschelBulkleyMu::castScalarModel
(
    const word& modelName
)
{
    if (modelName == muModel::typeName)
    {
        return dynamic_cast<muModel*>(this);
    }
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::HerschelBulkleyMu::muPatch
(
    const label patchi
) const
{
    const tmp<scalarField> sr(sMod_[sRateInd].boundaryField()[patchi]);
    return
    (
        min(mu0_, (tau0_ + k_*pow(sr(), n_))/(max(sr(), VSMALL)))
    );
}


Foam::tmp<Foam::scalarField> Foam::HerschelBulkleyMu::muInternal() const
{
    const tmp<scalarField> sr(sMod_[sRateInd].primitiveField());
    return
    (
        min(mu0_, (tau0_ + k_*pow(sr(), n_))/(max(sr(), VSMALL)))
    );
}


Foam::tmp<Foam::volScalarField>
Foam::HerschelBulkleyMu::muGeometric() const
{
    tmp<volScalarField> tMu
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh().time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("mu", muModel::modelDims, Zero)
        )
    );
    volScalarField& mu = tMu.ref();

    mu.primitiveFieldRef() = muInternal();
    forAll(mu.boundaryField(), patchi)
    {
        mu.boundaryFieldRef()[patchi].forceAssign(muPatch(patchi));
    }

    return tMu;
}


Foam::scalar Foam::HerschelBulkleyMu::muCell(const label celli) const
{
    const scalar sr(sMod_[sRateInd][celli]);
    return
    (
        min(mu0_, (tau0_ + k_*pow(sr, n_))/(max(sr, VSMALL)))
    );
}


bool Foam::HerschelBulkleyMu::read()
{
    k_ = dict_->lookup<scalar>("k");
    n_ = dict_->lookup<scalar>("n");
    tau0_ = dict_->lookup<scalar>("tau0");
    mu0_ = dict_->lookup<scalar>("mu0");

    return true;
}


// ************************************************************************* //
