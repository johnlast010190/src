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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/extrapolated/extrapolatedFvPatchScalarField.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

labelList extrapolatedFvPatchScalarField::fieldGradCalcTime_(0);
wordList extrapolatedFvPatchScalarField::fieldGradCalcName_(0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<vectorField> extrapolatedFvPatchScalarField::calcPatchCellToFace()
{
    tmp<vectorField> pF2C(new vectorField(this->patch().size(), vector::zero));

    const labelList& patchCells(this->patch().faceCells());
    const volVectorField& C(this->patch().boundaryMesh().mesh().C());

    forAll(patchCells, fI)
    {
        label cI = patchCells[fI];
        pF2C->operator[](fI) = C[cI];
    }

    pF2C.ref() = this->patch().Cf() - pF2C();

    return pF2C;
}

void extrapolatedFvPatchScalarField::calcPatchInternalGradient()
{
    word fname = this->internalField().name();
    word gradName = "extrapolatedFvPatchScalarField::grad" + fname;

    label fieldIndex = -1;
    forAll(fieldGradCalcName_, nI)
    {
        if (fieldGradCalcName_[nI] == fname)
        {
            fieldIndex = nI;
        }
    }

    if (fieldIndex == -1)
    {
        label cSize = fieldGradCalcName_.size();
        fieldGradCalcName_.setSize(cSize + 1);
        fieldGradCalcName_[cSize] = fname;
        fieldGradCalcTime_.setSize(cSize + 1);
        fieldGradCalcTime_[cSize] = -1;
        fieldIndex = cSize;
    }


    if
    (
        extrapolatedFvPatchScalarField::fieldGradCalcTime_[fieldIndex]
        != this->db().time().timeIndex()
    )
    {
        const volScalarField& f(this->db().lookupObject<volScalarField>(fname));

        autoPtr<volVectorField> gf
        (
            new volVectorField
            (
                IOobject
                (
                    gradName,
                    this->db().time().timeName(),
                    this->patch().boundaryMesh().mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::grad(f)
            )
        );

        gf->store(gf);

        extrapolatedFvPatchScalarField::fieldGradCalcTime_[fieldIndex]
            = this->db().time().timeIndex();
    }


    tmp<volVectorField> tfieldGrad
    (
        this->db().lookupObject<volVectorField>(gradName)
    );

    const labelList& patchCells(this->patch().faceCells());
    forAll(patchCells, fI)
    {
        label cI = patchCells[fI];
        patchInternalGradient_[fI] = tfieldGrad->operator[](cI);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedFvPatchScalarField::extrapolatedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    patchInternalGradient_(p.size(), vector::zero),
    patchCellToFace_(calcPatchCellToFace())
{}

extrapolatedFvPatchScalarField::extrapolatedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const Field<scalar>& f
)
:
    fvPatchField<scalar>(p, iF, f),
    patchInternalGradient_(p.size(), vector::zero),
    patchCellToFace_(calcPatchCellToFace())
{}

extrapolatedFvPatchScalarField::extrapolatedFvPatchScalarField
(
    const extrapolatedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper),
    patchInternalGradient_(ptf.patchInternalGradient_),
    patchCellToFace_(calcPatchCellToFace())
{}


extrapolatedFvPatchScalarField::extrapolatedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF, dict),
    patchInternalGradient_(p.size(), vector::zero),
    patchCellToFace_(calcPatchCellToFace())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(this->patchInternalField());
    }
}

extrapolatedFvPatchScalarField::extrapolatedFvPatchScalarField
(
    const extrapolatedFvPatchScalarField& ptf
)
:
    fvPatchField<scalar>(ptf),
    patchInternalGradient_(ptf.patchInternalGradient_),
    patchCellToFace_(calcPatchCellToFace())
{}

extrapolatedFvPatchScalarField::extrapolatedFvPatchScalarField
(
    const extrapolatedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(ptf, iF),
    patchInternalGradient_(ptf.patchInternalGradient_),
    patchCellToFace_(calcPatchCellToFace())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
    m(patchInternalGradient_, patchInternalGradient_);
}

void extrapolatedFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(ptf, addr);

    const extrapolatedFvPatchScalarField& fgptf =
        refCast<const extrapolatedFvPatchScalarField >(ptf);

    patchInternalGradient_.rmap(fgptf.patchInternalGradient_, addr);
}

void extrapolatedFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fvPatchScalarField::autoMapGIB(mapper);
    mapper.map(patchInternalGradient_, vector::zero);
}

// Update the coefficients associated with the patch field
void extrapolatedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->patch().boundaryMesh().mesh().changing())
    {
        patchCellToFace_ = calcPatchCellToFace();
    }

    calcPatchInternalGradient();

    fvPatchField<scalar>::updateCoeffs();
}

//- Return gradient at boundary
tmp<Field<scalar>> extrapolatedFvPatchScalarField::snGrad() const
{
    return (patchInternalGradient_ & this->patch().nf());
}


void extrapolatedFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const scalarField pif(this->patchInternalField());

    this->Field<scalar>::operator=
    (
        pif + (patchCellToFace_() & patchInternalGradient_)
    );

    fvPatchField<scalar>::evaluate();
}

tmp<Field<scalar>> extrapolatedFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>
        (new Field<scalar>(this->size(), pTraits<scalar>::one));
}


tmp<Field<scalar>> extrapolatedFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (patchCellToFace_() & patchInternalGradient_);
}


tmp<Field<scalar>> extrapolatedFvPatchScalarField::
gradientInternalCoeffs() const
{
    return tmp<Field<scalar>>
    (
        new Field<scalar>(this->size(), pTraits<scalar>::zero)
    );
}


tmp<Field<scalar>> extrapolatedFvPatchScalarField::
gradientBoundaryCoeffs() const
{
    return (patchInternalGradient_ & this->patch().nf());
}


// Write
void extrapolatedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, extrapolatedFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
