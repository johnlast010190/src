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
    (c) 2013 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "pointPatchFields/derived/constrainedOscillatingDisplacement/constrainedOscillatingDisplacementPointPatchVectorField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "global/constants/constants.H"
#include "fields/Fields/Field/SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constrainedOscillatingDisplacementPointPatchVectorField::
constrainedOscillatingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    amplitude_(0.0),
    frequency_(0.0),
    axis_(1,0,0),
    origin_(0,0,0),
    sourceRadius_(1),
    bufferWidth_(1),
    cubicF_()
{}


constrainedOscillatingDisplacementPointPatchVectorField::
constrainedOscillatingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    axis_(dict.lookupOrDefault<vector>
    (
        "axis",
        gSum(iF.mesh().GeoMesh<polyMesh>::operator()()
        .boundaryMesh()[p.index()].faceAreas())
        /gSum(mag(iF.mesh().GeoMesh<polyMesh>::operator()()
        .boundaryMesh()[p.index()].faceAreas()))
    )),
    origin_(dict.lookup("origin")),
    sourceRadius_(readScalar(dict.lookup("sourceRadius"))),
    bufferWidth_(readScalar(dict.lookup("bufferWidth"))),
    cubicF_()
{
    //calculate buffer function
    scalar x1 = sourceRadius_;
    scalar x12 = sqr(x1);
    scalar x13 = x12*x1;
    scalar x2 = sourceRadius_ + bufferWidth_;
    scalar x22 = sqr(x2);
    scalar x23 = x22*x2;

    scalar cf3Rcp = 3*x12*x2 - x23 -2*x13
        - 3.0/2.0 *(x12-x22)*(2*x1*x2-x22-x12)/(x1 - x2);

    cubicF_[3] = 1 / cf3Rcp;

    cubicF_[2] = -3.0/2.0*cubicF_[3]*(x12-x22)/(x1-x2);

    cubicF_[1] = -2*cubicF_[2]*x2 - 3*cubicF_[3]*x22;

    cubicF_[0] = 1 - cubicF_[1]*x1 - cubicF_[2]*x12 - cubicF_[3]*x13;

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


constrainedOscillatingDisplacementPointPatchVectorField::
constrainedOscillatingDisplacementPointPatchVectorField
(
    const constrainedOscillatingDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    sourceRadius_(ptf.sourceRadius_),
    bufferWidth_(ptf.bufferWidth_),
    cubicF_(ptf.cubicF_)
{}


constrainedOscillatingDisplacementPointPatchVectorField::
constrainedOscillatingDisplacementPointPatchVectorField
(
    const constrainedOscillatingDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    sourceRadius_(ptf.sourceRadius_),
    bufferWidth_(ptf.bufferWidth_),
    cubicF_(ptf.cubicF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void constrainedOscillatingDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    Field<vector>& pp(*this);

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    scalar angle = constant::mathematical::twoPi*frequency_*t.value();
    scalarField rp(mag(patch().localPoints() - origin_));
    scalar orad(sourceRadius_+bufferWidth_);

    pp.Field<vector>::operator=(axis_*amplitude_*sin(angle));

    forAll(rp, pI)
    {
        if (rp[pI] < sourceRadius_)
        {
            //do nothing
        }
        else if (rp[pI] > orad)
        {
            pp[pI] = vector::zero;
        }
        else
        {
            pp[pI] *= cubicF_.value(rp[pI]);

        }
    }


    fixedValuePointPatchField<vector>::updateCoeffs();
}


void constrainedOscillatingDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("amplitude", amplitude_);
    os.writeEntry("frequency", frequency_);

    os.writeEntry("axis", axis_);
    os.writeEntry("origin", origin_);

    os.writeEntry("sourceRadius", sourceRadius_);
    os.writeEntry("bufferWidth", bufferWidth_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    constrainedOscillatingDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
