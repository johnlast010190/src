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
    (c) 2016-2017 Wikki Ltd

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/derived/edgeNormalFixedValue/edgeNormalFixedValueFaPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/areaFields/areaFields.H"
#include "fields/faPatchFields/faPatchField/faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(p, iF),
    refValue_(p.size(), 0)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& ptf,
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchVectorField(ptf, p, iF, mapper),
    refValue_(mapper(ptf.refValue_))
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchVectorField(p, iF, dict),
    refValue_("refValue", dict, p.size())
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& pivpvf
)
:
    fixedValueFaPatchVectorField(pivpvf),
    refValue_(pivpvf.refValue_)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& pivpvf,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(pivpvf, iF),
    refValue_(pivpvf.refValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::edgeNormalFixedValueFaPatchVectorField::autoMap
(
    const faPatchFieldMapper& m
)
{
    fixedValueFaPatchVectorField::autoMap(m);
    m(refValue_, refValue_);
}


void Foam::edgeNormalFixedValueFaPatchVectorField::rmap
(
    const faPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFaPatchVectorField::rmap(ptf, addr);

    const edgeNormalFixedValueFaPatchVectorField& tiptf =
        refCast<const edgeNormalFixedValueFaPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::edgeNormalFixedValueFaPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    forceAssign(refValue_*patch().edgeNormals());
}


void Foam::edgeNormalFixedValueFaPatchVectorField::write(Ostream& os) const
{
    fixedValueFaPatchVectorField::write(os);
    refValue_.writeEntry("refValue", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeFaPatchTypeField
(
    faPatchVectorField,
    edgeNormalFixedValueFaPatchVectorField
);

} // End namespace Foam


// ************************************************************************* //
