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
    (c) 2011-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/unidirectionalVelocityOutlet/unidirectionalVelocityOutletFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unidirectionalVelocityOutletFvPatchVectorField::unidirectionalVelocityOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(p, iF),
    sink_(1e+9)
{}


unidirectionalVelocityOutletFvPatchVectorField::unidirectionalVelocityOutletFvPatchVectorField
(
    const unidirectionalVelocityOutletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchVectorField(ptf, p, iF, mapper),
    sink_(ptf.sink_)
{}


unidirectionalVelocityOutletFvPatchVectorField::unidirectionalVelocityOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchVectorField(p, iF),
    sink_(dict.lookupOrDefault<scalar>("Cs", 1e+9))
{
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    fvPatchVectorField::forceAssign(vectorField("value", dict, p.size()));
}


unidirectionalVelocityOutletFvPatchVectorField::unidirectionalVelocityOutletFvPatchVectorField
(
    const unidirectionalVelocityOutletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(ptf, iF),
    sink_(ptf.sink_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void unidirectionalVelocityOutletFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Field<scalar>& phip =
        this->patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
        );

    scalarField upn(*this & patch().nf());
    this->valueFraction() = 1.0 - min(pos0(upn), pos0(phip));

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::unidirectionalVelocityOutletFvPatchVectorField::manipulateMatrix
(
    fvVectorMatrix& matrix
)
{
/*

    //Add explicit sink to momentum equation to prevent inflow
    // in cells adjacent to patch
    const surfaceScalarField phi(db().lookupObject<surfaceScalarField>(phiName_));
    label index = patch().index();
    const fvsPatchField<scalar>& phip(phi.boundaryField()[index]);
    const fvsPatchField<scalar>& phipOld
    (
        phi.oldTime().boundaryField()[index]
    );
    scalarField upn = (*this & patch().nf());


    const labelList& faceCells = this->patch().faceCells();
    const scalarField& V = patch().boundaryMesh().mesh().V();

    scalarField& UaDiag = matrix.diag();

    forAll(faceCells, fcI)
    {
        label cI = faceCells[fcI];

        UaDiag[cI] += V[cI] * sink_ * max(neg(phip[fcI]), neg(upn[fcI]));
        UaDiag[cI] += V[cI] * 1.0 * neg(phipOld[fcI]);
    }

    const surfaceScalarField phi(db().lookupObject<surfaceScalarField>(phiName_));
    label index = patch().index();
    const fvsPatchField<scalar>& phip(phi.boundaryField()[index]);
    scalarField upn = (*this & patch().nf());

    const labelList& faceCells = this->patch().faceCells();
    labelList reversingCells(faceCells.size(), -1);
    label nBlockedFaces = 0;

    forAll(faceCells, fcI)
    {
        if (phip[fcI] <= 0 || upn[fcI] <= 0)
        {
            reversingCells[nBlockedFaces] = faceCells[fcI];
            nBlockedFaces++;
        }
    }
    reversingCells.setSize(nBlockedFaces);
    matrix.setValues(reversingCells, vectorField(nBlockedFaces, vector::zero));
*/
}


// Write
void unidirectionalVelocityOutletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<scalar>(os, "Cs", 1e9, sink_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchVectorField,
    unidirectionalVelocityOutletFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
