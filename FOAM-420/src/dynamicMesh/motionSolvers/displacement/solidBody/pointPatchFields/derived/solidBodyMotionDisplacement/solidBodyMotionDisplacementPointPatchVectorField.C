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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolvers/displacement/solidBody/pointPatchFields/derived/solidBodyMotionDisplacement/solidBodyMotionDisplacementPointPatchVectorField.H"
#include "fields/Fields/transformField/transformField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    isSolidBody_(false),
    SBMFPtr_(),
    coorFramePtr_(nullptr)
{}


solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict, false),
    isSolidBody_(dict.found("solidBodyMotionFunction")),
    SBMFPtr_
    (
        isSolidBody_
      ? solidBodyMotionFunction::New(dict, this->db().time())
      : nullptr
    ),
    coorFramePtr_(nullptr)
{
    if (!isSolidBody_)
    {
        coorFramePtr_ =
                coordinateFrame::lookupNew
                (
                    dynamic_cast<const fvMesh&>(this->db()),
                    dict
                );

        if (!coorFramePtr_->anyDynamic())
        {
            FatalErrorInFunction
                << "coordinateFrame " << coorFramePtr_->name()
                << " is linked to dynamic mesh but it is not dynamic."
                << abort(FatalError);
        }
    }
    if (!dict.found("value"))
    {
        // Necessary to support old solvers
        // (foamSolve has it's own updateStates)
        coordinateFrame::updateStates(patch().boundaryMesh().mesh().thisDb());
        const septernion transform =
            isSolidBody_
          ? SBMFPtr_().transformation()
          : coorFramePtr_->transformation();

        // Determine current local points and offset
        fixedValuePointPatchVectorField::forceAssign
        (
            transformPoints(transform, localPoints0()) - localPoints0()
        );
    }
}


solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const solidBodyMotionDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    isSolidBody_(ptf.isSolidBody_),
    SBMFPtr_(isSolidBody_ ? ptf.SBMFPtr_().clone().ptr() : nullptr),
    coorFramePtr_(!isSolidBody_ ? ptf.coorFramePtr_ : nullptr)
{
    // Necessary to support old solvers
    // (foamSolve has it's own updateStates)
    coordinateFrame::updateStates(patch().boundaryMesh().mesh().thisDb());

    // For safety re-evaluate
    const septernion transform =
        isSolidBody_
      ? SBMFPtr_().transformation()
      : coorFramePtr_->transformation();

    fixedValuePointPatchVectorField::forceAssign
    (
        transformPoints(transform, localPoints0()) - localPoints0()
    );
}


solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const solidBodyMotionDisplacementPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    isSolidBody_(ptf.isSolidBody_),
    SBMFPtr_(isSolidBody_ ? ptf.SBMFPtr_().clone().ptr() : nullptr),
    coorFramePtr_(!isSolidBody_ ? ptf.coorFramePtr_ : nullptr)
{}


solidBodyMotionDisplacementPointPatchVectorField::
solidBodyMotionDisplacementPointPatchVectorField
(
    const solidBodyMotionDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    isSolidBody_(ptf.isSolidBody_),
    SBMFPtr_(isSolidBody_ ? ptf.SBMFPtr_().clone().ptr() : nullptr),
    coorFramePtr_(!isSolidBody_ ? ptf.coorFramePtr_ : nullptr)
{
    // Necessary to support old solvers (foamSolve has it's own updateStates)
    coordinateFrame::updateStates(patch().boundaryMesh().mesh().thisDb());

    // For safety re-evaluate
    const septernion transform =
        isSolidBody_
      ? SBMFPtr_().transformation()
      : coorFramePtr_->transformation();

    fixedValuePointPatchVectorField::forceAssign
    (
        transformPoints(transform, localPoints0()) - localPoints0()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const pointField&
solidBodyMotionDisplacementPointPatchVectorField::localPoints0() const
{
    if (!localPoints0Ptr_.valid())
    {
        pointIOField points0
        (
            IOobject
            (
                "points",
                this->db().time().constant(),
                polyMesh::meshSubDir,
                this->db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        localPoints0Ptr_.reset(new pointField(points0, patch().meshPoints()));
    }
    return localPoints0Ptr_();
}


void solidBodyMotionDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Necessary to support old solvers
    // (foamSolve has it's own updateStates)
    coordinateFrame::updateStates(patch().boundaryMesh().mesh().thisDb());

    const septernion transform =
        isSolidBody_
      ? SBMFPtr_().transformation()
      : coorFramePtr_->transformation();

    // Determine current local points and offset
    fixedValuePointPatchVectorField::forceAssign
    (
        transformPoints(transform, localPoints0())
       -localPoints0()
    );

    fixedValuePointPatchVectorField::updateCoeffs();
}


void solidBodyMotionDisplacementPointPatchVectorField::write(Ostream& os) const
{
    // Note: write value
    fixedValuePointPatchVectorField::write(os);

    if (isSolidBody_)
    {
        os.writeEntry
        (
            solidBodyMotionFunction::typeName,
            SBMFPtr_->type()
        );
        os  << indent << word(SBMFPtr_->type() + "Coeffs");
        SBMFPtr_->writeData(os);
    }
    else
    {
        os.writeEntry("referenceFrame", coorFramePtr_->name());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    solidBodyMotionDisplacementPointPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
