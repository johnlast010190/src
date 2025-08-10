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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fixedMultiPhaseHeatFluxFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "twoPhaseSystem.H"
#include "PhaseSystems/ThermalPhaseChangePhaseSystem/ThermalPhaseChangePhaseSystem.H"
#include "PhaseSystems/MomentumTransferPhaseSystem/MomentumTransferPhaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity/ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel/PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    relax_(1.0)
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    q_("q", dict, p.size()),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0))
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fixedMultiPhaseHeatFluxFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    q_(mapper(psf.q_)),
    relax_(psf.relax_)
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fixedMultiPhaseHeatFluxFvPatchScalarField& psf
)
:
    fixedValueFvPatchScalarField(psf),
    q_(psf.q_),
    relax_(psf.relax_)
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fixedMultiPhaseHeatFluxFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    q_(psf.q_),
    relax_(psf.relax_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the fluid model
    const ThermalPhaseChangePhaseSystem
    <
        MomentumTransferPhaseSystem<twoPhaseSystem>
    >& fluid =
        refCast
        <
            const ThermalPhaseChangePhaseSystem
            <
                MomentumTransferPhaseSystem<twoPhaseSystem>
            >
        >
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const scalarField& Tp = *this;

    scalarField A(Tp.size(), scalar(0));
    scalarField B(Tp.size(), scalar(0));
    scalarField Q(Tp.size(), scalar(0));

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        const fluidThermo& thermo = phase.thermo();

        const fvPatchScalarField& alpha =
            phase.boundaryField()[patch().index()];

        const fvPatchScalarField& T =
            thermo.T().boundaryField()[patch().index()];

        const scalarField kappaEff
        (
            thermo.kappaEff
            (
                phase.turbulence().alphat(patch().index()),
                patch().index()
            )
        );

        if (debug)
        {
            scalarField q0(T.snGrad()*alpha*kappaEff);
            Q += q0;

            Info<< patch().name() << " " << phase.name()
                << ": Heat flux " << gMin(q0) << " - " << gMax(q0) << endl;
        }

        A += T.patchInternalField()*alpha*kappaEff*patch().deltaCoeffs();
        B += alpha*kappaEff*patch().deltaCoeffs();
    }

    if (debug)
    {
        Info<< patch().name() << " " << ": overall heat flux "
            << gMin(Q) << " - " << gMax(Q) << endl;
    }

    forceAssign((1 - relax_)*Tp + relax_*(q_ + A)/(B));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("relax", relax_);
    q_.writeEntry("q", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedMultiPhaseHeatFluxFvPatchScalarField
    );
}


// ************************************************************************* //
