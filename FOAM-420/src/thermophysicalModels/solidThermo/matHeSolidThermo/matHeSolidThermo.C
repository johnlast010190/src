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
    (c) 2022 Esi Ltd.
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matHeSolidThermo.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::matHeSolidThermo<BasicSolidThermo, MixtureType>::calculate()
{
    matScalarTable& sModels =
        this->materials_.sTable
        (
            this->phaseName_,
            word::null
        );
    const baseModels<scalar>& TMod = (*sModels[TModel::typeName]);
    const baseModels<scalar>& rhoMod = (*sModels[rhoModel::typeName]);
    const baseModels<scalar>& kappaMod = (*sModels[kappaModel::typeName]);
    const baseModels<scalar>& CpvMod = (*sModels[CpvModel::typeName]);
    const baseModels<scalar>& HEMod = (*sModels[HEModel::typeName]);

    this->T_.primitiveFieldRef() = TMod.primitiveField();
    this->rho_.primitiveFieldRef() = rhoMod.primitiveField();
    this->alpha_.primitiveFieldRef() = kappaMod.primitiveField()/CpvMod.primitiveField();

    forAll(this->T_.boundaryField(), patchi)
    {
        if (this->T_.boundaryFieldRef()[patchi].fixesValue())
        {
            this->he().boundaryFieldRef()[patchi].forceAssign(
                HEMod.boundaryField()[patchi]
            );
        }
        else
        {
            this->T_.boundaryFieldRef()[patchi].forceAssign(
                TMod.boundaryField()[patchi]
            );
        }
        this->rho_.boundaryFieldRef()[patchi].forceAssign(
            rhoMod.boundaryField()[patchi]
        );
        this->alpha_.boundaryFieldRef()[patchi].forceAssign(
            kappaMod.boundaryField()[patchi]/CpvMod.boundaryField()[patchi]
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::matHeSolidThermo<BasicSolidThermo, MixtureType>::
matHeSolidThermo
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    matHeThermo<BasicSolidThermo, MixtureType>(obr, phaseName)
{
    calculate();
}


template<class BasicSolidThermo, class MixtureType>
Foam::matHeSolidThermo<BasicSolidThermo, MixtureType>::
matHeSolidThermo
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
:
    matHeThermo<BasicSolidThermo, MixtureType>(obr, dict, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::matHeSolidThermo<BasicSolidThermo, MixtureType>::~matHeSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::matHeSolidThermo<BasicSolidThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volVectorField>
Foam::matHeSolidThermo<BasicSolidThermo, MixtureType>::Kappa() const
{
    return
        this->materials_.vTable
        (
            this->phaseName_,
            word::null
        )[vKappaModel::typeName]->operator()();
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::vectorField>
Foam::matHeSolidThermo<BasicSolidThermo, MixtureType>::Kappa
(
    const label patchi
) const
{
    return
        this->materials_.vTable
        (
            this->phaseName_,
            word::null
        )[vKappaModel::typeName]->boundaryField()[patchi];
}


// ************************************************************************* //
