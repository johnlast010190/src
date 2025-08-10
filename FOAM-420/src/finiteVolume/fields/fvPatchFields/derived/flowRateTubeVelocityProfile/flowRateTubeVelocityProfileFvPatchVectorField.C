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
    (c) 2016 OpenFOAM Foundation
    (c) 2017 Vascular Flow
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/flowRateTubeVelocityProfile/flowRateTubeVelocityProfileFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateTubeVelocityProfileFvPatchVectorField::
flowRateTubeVelocityProfileFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowRateBase(patch(), db(), false),
    n_(2.0)
{}


Foam::flowRateTubeVelocityProfileFvPatchVectorField::
flowRateTubeVelocityProfileFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowRateBase(patch(), db(), dict, false),
    n_(dict.lookupOrDefault<scalar>("n", 2.0))
{
    if (n_ <= 0.0)
    {
        FatalIOErrorIn
        (
            "flowRateTubeVelocityProfileFvPatchVectorField::"
            "flowRateTubeVelocityProfileFvPatchVectorField"
            "(const fvPatch&, const DimensionedField<vector, volMesh>&,"
            " const dictionary&)",
            dict
        )   << "The degree of polynomial 'n' must be positive"
            << exit(FatalIOError);
    }

    // Value field require if mass based
    if (dict.found("value"))
    {
        // Initialise with the value entry
        forceAssign(vectorField("value", dict, p.size()));
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::flowRateTubeVelocityProfileFvPatchVectorField::
flowRateTubeVelocityProfileFvPatchVectorField
(
    const flowRateTubeVelocityProfileFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flowRateBase(patch(), db(), ptf, false),
    n_(ptf.n_)
{}


Foam::flowRateTubeVelocityProfileFvPatchVectorField::
flowRateTubeVelocityProfileFvPatchVectorField
(
    const flowRateTubeVelocityProfileFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    flowRateBase(patch(), db(), ptf, false),
    n_(ptf.n_)
{}


Foam::flowRateTubeVelocityProfileFvPatchVectorField::
flowRateTubeVelocityProfileFvPatchVectorField
(
    const flowRateTubeVelocityProfileFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    flowRateBase(patch(), db(), ptf, false),
    n_(ptf.n_)
{
    // Evaluate the profile if defined.
    // (not sure what this does, it was copied from
    //  uniformFixedValueFvPatchField.C)
    //if (ptf.flowRate_.valid())
    //{
        // this call causes issues with caseSetup potential flow
        // initialisation as rho is not yet defined while reading U!
    //    fixedValueFvPatchVectorField::evaluate();
    //}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::flowRateTubeVelocityProfileFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<vectorField> vic(new vectorField(this->size(), pTraits<vector>::one));
    vic.ref() *= neg(currentFlowRate());

    return vic;
}


Foam::tmp<Foam::vectorField>
Foam::flowRateTubeVelocityProfileFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<vectorField> vbc(new vectorField(*this));
    vbc.ref() *= pos0(currentFlowRate());

    return vbc;
}


void Foam::flowRateTubeVelocityProfileFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Circular patch area
    scalar A = gSum(patch().magSf());

    // Circular patch radius is the weighted mean of patch face centres
    vector O = gSum(patch().magSf()*patch().Cf())/A;

    // Circular patch normal vector is the weighted mean of patch face normal
    // vectors
    vector n = gSum(patch().magSf()*patch().nf())/A;

    // The patch must be a circular patch for this to be true!
    scalar R = sqrt(A/constant::mathematical::pi);

    // Radial position of cell centres
    tmp<scalarField> r = mag(patch().Cf() - O);

    tmp<vectorField> tUp
    (
       -n*currentFlowRate()/A
       *(n_ + 2.0)/n_*(1.0 - pow((r/R), n_))
    );
    vectorField& Up = tUp.ref();

    Up /= multiplyByRhoOrOne(scalarField(this->size(), 1));

    frameFieldUpdate(Up);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateTubeVelocityProfileFvPatchVectorField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchVectorField::write(os);
    flowRateBase::write(os);
    os.writeEntry("n", n_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        flowRateTubeVelocityProfileFvPatchVectorField
    );
}

// ************************************************************************* //
