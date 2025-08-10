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
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "CassonMu.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CassonMu, 0);
    addToRunTimeSelectionTable(materialModel, CassonMu, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CassonMu::CassonMu
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
    CassonMu::read();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::CassonMu::updateTable(const word& modelName)
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


Foam::baseModels<Foam::scalar>* Foam::CassonMu::castScalarModel
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


Foam::tmp<Foam::scalarField> Foam::CassonMu::muPatch
(
    const label patchi
) const
{
    return
        max
        (
            muMin_,
            min
            (
                muMax_,
                sqr
                (
                    sqrt
                    (
                        tau0_
                       /max
                        (
                            sMod_[sRateInd].boundaryField()[patchi],
                            VSMALL
                        )
                    ) + sqrt(m_)
                )
            )
        );
}


Foam::tmp<Foam::scalarField> Foam::CassonMu::muInternal() const
{
    return
        max
        (
            muMin_,
            min
            (
                muMax_,
                sqr
                (
                    sqrt
                    (
                        tau0_
                       /max
                        (
                            sMod_[sRateInd].primitiveField(),
                            VSMALL
                        )
                    ) + sqrt(m_)
                )
            )
        );
}


Foam::tmp<Foam::volScalarField> Foam::CassonMu::muGeometric() const
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


Foam::scalar Foam::CassonMu::muCell(const label celli) const
{
    return
        max
        (
            muMin_,
            min
            (
                muMax_,
                sqr
                (
                    sqrt
                    (
                        tau0_
                       /max(sMod_[sRateInd][celli], VSMALL)
                    ) + sqrt(m_)
                )
            )
        );
}


bool Foam::CassonMu::read()
{
    m_ = dict_->lookup<scalar>("m");
    tau0_ = dict_->lookup<scalar>("tau0");
    muMin_ = dict_->lookup<scalar>("muMin");
    muMax_ = dict_->lookup<scalar>("muMax");

    return true;
}


// ************************************************************************* //
