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
    (c) 2010-2012 Esi Ltd.
    (c) 2006 Mark Olesen

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/wheelVelocity/wheelVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void wheelVelocityFvPatchVectorField::calcVelocity()
{
    if (size())
    {
        if (mag(axis_) == 0)
        {
            FatalErrorInFunction
                << "Zero sized axis specified for wheelVelocity." << endl;
        }
        vector naxis = axis_/mag(axis_);

        if (mag(contactRadius_) == 0)
        {
            FatalErrorInFunction
                << "Zero contact radius specified for wheelVelocity." << endl;
        }

        vector Omega = naxis * hubSpeed_/contactRadius_;

        vectorField origToCf (patch().Cf() - origin_);

        vectorField& Up = *this;

        Up = (Omega ^ origToCf);

        //no surface normal motion
        vectorField nf( patch().nf() );
        Up -= nf* (Up & nf);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wheelVelocityFvPatchVectorField::
wheelVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(vector::zero),
    axis_(vector(1,0,0)),
    hubSpeed_(0),
    contactRadius_(1)
{}

wheelVelocityFvPatchVectorField::
wheelVelocityFvPatchVectorField
(
    const wheelVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    hubSpeed_(ptf.hubSpeed_),
    contactRadius_(ptf.contactRadius_)
{
    calcVelocity();
}

wheelVelocityFvPatchVectorField::
wheelVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    hubSpeed_(readScalar(dict.lookup("hubSpeed"))),
    contactRadius_(readScalar(dict.lookup("contactRadius")))
{
    calcVelocity();
}

wheelVelocityFvPatchVectorField::
wheelVelocityFvPatchVectorField
(
    const wheelVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    hubSpeed_(ptf.hubSpeed_),
    contactRadius_(ptf.contactRadius_)
{
    calcVelocity();
}

wheelVelocityFvPatchVectorField::
wheelVelocityFvPatchVectorField
(
    const wheelVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    hubSpeed_(ptf.hubSpeed_),
    contactRadius_(ptf.contactRadius_)
{
    calcVelocity();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wheelVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void wheelVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    os.writeEntry("hubSpeed", hubSpeed_);
    os.writeEntry("contactRadius", contactRadius_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    wheelVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
