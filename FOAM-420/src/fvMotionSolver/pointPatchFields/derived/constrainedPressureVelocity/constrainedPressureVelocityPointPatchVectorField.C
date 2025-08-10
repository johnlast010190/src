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

#include "pointPatchFields/derived/constrainedPressureVelocity/constrainedPressureVelocityPointPatchVectorField.H"
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

scalar constrainedPressureVelocityPointPatchVectorField::lineGrad
(
    const Function1<scalar>& f,
    scalar x,
    scalar dx
)
{
    return ((f.value(x+dx) - f.value(x-dx))/(2*dx));
}

constrainedPressureVelocityPointPatchVectorField::
constrainedPressureVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    pressure_(),
    rho0_(1.0),
    gamma_(0.0),
    p0_(0.0),
    axis_(1,0,0),
    origin_(0,0,0),
    sourceRadius_(1),
    bufferWidth_(1),
    cubicF_()
{}


constrainedPressureVelocityPointPatchVectorField::
constrainedPressureVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    pressure_(Function1<scalar>::New("pressure", dict)),
    rho0_(readScalar(dict.lookup("rho0"))),
    gamma_(readScalar(dict.lookup("gamma"))),
    p0_(readScalar(dict.lookup("p0"))),
    axis_(vector::zero),
    origin_(dict.lookup("origin")),
    sourceRadius_(readScalar(dict.lookup("sourceRadius"))),
    bufferWidth_(readScalar(dict.lookup("bufferWidth"))),
    cubicF_()
{
    if (bufferWidth_/sourceRadius_ > SMALL)
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

    }
    else
    {
        cubicF_[3] = 0;

        cubicF_[2] = 0;

        cubicF_[1] = 0;

        cubicF_[0] = 0;
    }

    //calculate the mean normal vector in the loudspeaker zone
    const polyPatch& polp = iF.mesh().GeoMesh<polyMesh>::operator()()
            .boundaryMesh()[p.index()];

    const vectorField::subField& faceCentres(polp.faceCentres());
    const vectorField::subField& faceAreas(polp.faceAreas());
    scalarField rp(mag(faceCentres - origin_));

    forAll(rp, i)
    {
        if (rp[i] < sourceRadius_)
        {
            axis_ += faceAreas[i];
        }
    }

    reduce(axis_, sumOp<vector>());

    axis_ /= mag(axis_);

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


constrainedPressureVelocityPointPatchVectorField::
constrainedPressureVelocityPointPatchVectorField
(
    const constrainedPressureVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    pressure_
    (
        ptf.pressure_.valid()
      ? ptf.pressure_().clone().ptr()
      : nullptr
    ),
    rho0_(ptf.rho0_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    sourceRadius_(ptf.sourceRadius_),
    bufferWidth_(ptf.bufferWidth_),
    cubicF_(ptf.cubicF_)
{}


constrainedPressureVelocityPointPatchVectorField::
constrainedPressureVelocityPointPatchVectorField
(
    const constrainedPressureVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    pressure_
    (
        ptf.pressure_.valid()
      ? ptf.pressure_().clone().ptr()
      : nullptr
    ),
    rho0_(ptf.rho0_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    sourceRadius_(ptf.sourceRadius_),
    bufferWidth_(ptf.bufferWidth_),
    cubicF_(ptf.cubicF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void constrainedPressureVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();

    scalarField rp(mag(patch().localPoints() - origin_));
    scalar orad(sourceRadius_+bufferWidth_);

    Field<vector>& Up(*this);
    scalar c0 = sqrt(gamma_*p0_/rho0_);

    //if (this->db().time().timeIndex() == 1)
    {
        //straight up
        Up = -axis_*1/(rho0_*c0)*(pressure_->value(t));
    }
    /*else //commented out since change was negligible
    {
        //RK4
        scalar dt = this->db().time().deltaT().value();

        scalar k1 = lineGrad(pressure_(), t-dt, 0.1*dt);

        scalar k2 = lineGrad(pressure_(), t-0.5*dt, 0.1*dt);

        scalar k3 = lineGrad(pressure_(), t, 0.1*dt);

        Up -= axis_ * 1/(rho0_*c0) * (1.0/6.0) * dt * (k1 + 4*k2 + k3);

        Info<< "Current/old method Up_mag: " << (Up[0] & axis_  )
             << " " << -1/(rho0_*c0)*(pressure_->value(t)) << endl;

    }*/


    forAll(rp, pI)
    {
        if (rp[pI] < sourceRadius_)
        {
            //do nothing
        }
        else if (rp[pI] > orad)
        {
            Up[pI] = vector::zero;
        }
        else
        {
            Up[pI] *= cubicF_.value(rp[pI]);

        }
    }

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void constrainedPressureVelocityPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("rho0", rho0_);
    os.writeEntry("gamma", gamma_);
    os.writeEntry("p0", p0_);

    os.writeEntry("origin", origin_);

    os.writeEntry("sourceRadius", sourceRadius_);
    os.writeEntry("bufferWidth", bufferWidth_);

    pressure_->writeData(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    constrainedPressureVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
