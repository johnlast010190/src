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
    (c) 2016-2021 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/


#include "BirdCarreauMu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BirdCarreauMu, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        BirdCarreauMu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BirdCarreauMu::BirdCarreauMu
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
    BirdCarreauMu::read();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::BirdCarreauMu::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>* Foam::BirdCarreauMu::castScalarModel
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


Foam::tmp<Foam::scalarField> Foam::BirdCarreauMu::muPatch
(
    const label patchi
) const
{
    return
        muInf_
      + (mu0_ - muInf_)
       *pow
        (
            1.0 + pow(k_*sMod_[sRateInd].boundaryField()[patchi], a_),
            (n_ - 1.0)/a_
        );
}


Foam::tmp<Foam::scalarField> Foam::BirdCarreauMu::muInternal() const
{
    return
        muInf_
      + (mu0_ - muInf_)
       *pow
        (
            1.0 + pow(k_*sMod_[sRateInd].primitiveField(), a_),
            (n_ - 1.0)/a_
        );
}


Foam::tmp<Foam::volScalarField>
Foam::BirdCarreauMu::muGeometric() const
{
    dimensionedScalar muInf("muInf", dimDynamicViscosity, muInf_);
    dimensionedScalar mu0("mu0", dimDynamicViscosity, mu0_);
    dimensionedScalar a("a", dimless, a_);
    dimensionedScalar n("n", dimless, n_);
    dimensionedScalar k("k", dimTime, k_);

    return
        muInf
      + (mu0 - muInf)
       *pow(1.0 + pow(k*sMod_[sRateInd](), a),(n - 1.0)/a);
}


Foam::scalar Foam::BirdCarreauMu::muCell(const label celli) const
{
    return
        muInf_
      + (mu0_ - muInf_)
       *pow
        (
            1.0 + pow(k_*sMod_[sRateInd][celli], a_),
            (n_ - 1.0)/a_
        );
}


bool Foam::BirdCarreauMu::read()
{
    mu0_ = dict_->lookup<scalar>("mu0");
    muInf_ = dict_->lookup<scalar>("muInf");
    k_ = dict_->lookup<scalar>("k");
    n_ = dict_->lookup<scalar>("n");
    a_ = dict_->lookupOrDefault<scalar>("a", 2.0);

    return true;
}


// ************************************************************************* //
