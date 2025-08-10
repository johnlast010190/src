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

\*---------------------------------------------------------------------------*/

#include "solidChemistryModel/solidChemistryModel.H"
#include "mixtures/reactingMixture/reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::solidChemistryModel<CompType, SolidThermo>::
solidChemistryModel
(
    const objectRegistry& obr,
    const word& phaseName
)
:
    CompType(obr, phaseName),
    ODESystem(),
    Ys_(this->solidThermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<SolidThermo>& >
        (
            this->solidThermo()
        )
    ),
    solidThermo_
    (
        dynamic_cast<const reactingMixture<SolidThermo>& >
        (
            this->solidThermo()
        ).speciesData()
    ),
    nSolids_(Ys_.size()),
    nReaction_(reactions_.size()),
    RRs_(nSolids_),
    reactingCells_(this->mesh().nCells(), true)
{
    // create the fields for the chemistry sources
    forAll(RRs_, fieldi)
    {
        RRs_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRs." + Ys_[fieldi].name(),
                    obr.time().timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
   }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::solidChemistryModel<CompType, SolidThermo>::
~solidChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::scalar Foam::solidChemistryModel<CompType, SolidThermo>::solve
(
    const scalarField& deltaT
)
{
    NotImplemented;
}


template<class CompType, class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<CompType, SolidThermo>::tc() const
{
    NotImplemented;
}


template<class CompType, class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<CompType, SolidThermo>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Ys_, i)
        {
            forAll(Qdot, celli)
            {
                scalar hf = solidThermo_[i].Hc();
                Qdot[celli] -= hf*RRs_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class CompType, class SolidThermo>
void Foam::solidChemistryModel<CompType, SolidThermo>::setCellReacting
(
    const label celli,
    const bool active
)
{
    reactingCells_[celli] = active;
}

// ************************************************************************* //
