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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "matHePsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::matHePsiThermo<BasicPsiThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        this->materials_.updateScalarField(this->phasePropertyName("p"), p.oldTime());
        this->materials_.updateScalarField(this->phasePropertyName("T"), T.oldTime());
        this->materials_.updateScalarField(this->phasePropertyName(he.name()), he.oldTime());
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    // Update temperature internal field
    if (this->solvesTFromhe())
    {
        T.primitiveFieldRef() =
            this->materials_
            (
                TModel::typeName, this->phaseName_
            ).primitiveField();
    }
    else
    {
        he.primitiveFieldRef() =
            this->materials_
            (
                HEModel::typeName, this->phaseName_
            ).primitiveField();
    }

    // Is temperature constant?
    const bool isTConst =
        obr_.found("TRef")
      ? obr_.lookupObject<refScalarField>("TRef").isConst()
      : false;

    // Update he and T boundary fields
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    forAll(TBf, patchi)
    {
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];

        if (this->solvesTFromhe())
        {
            if (pT.fixesValue() && !isTConst)
            {
                phe.forceAssign(
                    this->materials_
                    (
                        HEModel::typeName,
                        this->phaseName_
                    ).boundaryField()[patchi]
                );
            }
            else
            {
                pT.forceAssign(
                    this->materials_
                    (
                        TModel::typeName,
                        this->phaseName_
                    ).boundaryField()[patchi]
                );
            }
        }
        else
        {
            phe.forceAssign(
                this->materials_
                (
                    HEModel::typeName,
                    this->phaseName_
                ).boundaryField()[patchi]
            );
        }
    }
    psi.forceAssign(this->materials_(psiModel::typeName, this->phaseName_)());
    mu.forceAssign(this->materials_(muModel::typeName, this->phaseName_)());
    alpha.forceAssign(this->materials_(alphahModel::typeName, this->phaseName_)());

    // Update back when old times is true and we are doing current timestep
    // which isn't 0 (nOldTimes() = 0 means there were no old times and update isn't necessary)
    if (doOldTimes && this->p_.nOldTimes() && (this->p_.nOldTimes() == p.nOldTimes()))
    {
        this->materials_.updateScalarField(this->phasePropertyName("p"), this->p_);
        this->materials_.updateScalarField(this->phasePropertyName("T"), this->T_);
        this->materials_.updateScalarField(this->phasePropertyName(he.name()), this->he_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::matHePsiThermo<BasicPsiThermo, MixtureType>::matHePsiThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    matHeThermo<BasicPsiThermo, MixtureType>(obr, phaseName),
    obr_(obr)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::matHePsiThermo<BasicPsiThermo, MixtureType>::~matHePsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::matHePsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    BasicPsiThermo::correct();
}


// ************************************************************************* //
