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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/motionCoordinateFrame/motionCoordinateFrame.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionCoordinateFrame, 0);
    addToRunTimeSelectionTable
    (
        coordinateFrame,
        motionCoordinateFrame,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::motionCoordinateFrame::dimensionCheck
(
    const dimensionSet& dim1,
    const dimensionSet& dim2
) const
{
    if (dim1 != dim2)
    {
        FatalErrorInFunction
            << "UEqn dimensions do not match " << dim2
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionCoordinateFrame::motionCoordinateFrame
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& frameName
)
:
    coordinateFrame(mesh, dict, frameName),
    motionFunction_(nullptr)
{
    const dictionary& motionDict = dict.subDict("motionFunction");
    motionFunction_ =
        solidBodyMotionFunction::New(mesh_, motionDict, frameName);

    // Should be defined from the motion so sub-dicts work correctly
    this->isIncrementalMotion() = motionFunction_().isIncrementalMotion();

    // Check if motion allows outer corrector motion
    this->outerCorrectorMotion() = motionFunction_().outerCorrectorMotion();

    // Initialize the coordinate frame state data
    motionCoordinateFrame::updateState();
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::motionCoordinateFrame::~motionCoordinateFrame()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::motionCoordinateFrame::updateState() const
{
    if (!isUpdated())
    {
        // First store old time data, then update
        storeOldTimes();

        // Update current coordinate system
        const label index = nFrameCorrector();
        const label size = index + 1;
        transformations_.setSize(size);
        if (!outerCorrectorMotion() && index > 0)
        {
            transformations_[index] = septernion::I;
        }
        else
        {
            transformations_[index] = motionFunction_().transformation();
        }

        updateCoordinateSystem();

        // Update the coordinate frame state data
        velocities_.setSize(size);
        if (!outerCorrectorMotion() && index > 0)
        {
            velocities_[index] = vectorTuple(Zero, Zero);
        }
        else
        {
            velocities_[index] = motionFunction_().velocity();
        }

        if (mesh_.thisDb().time().timeIndex() == 0)
        {
            const_cast<vectorTuple&>(oldTime().velocity()) =
                velocities_[index];
        }

        accelerations_.setSize(size);
        if (!outerCorrectorMotion() && index > 0)
        {
            accelerations_[index] = vectorTuple(Zero, Zero);
        }
        else
        {
            accelerations_[index] = motionFunction_().acceleration();
        }

        updateIndex_ = obr_.time().timeIndex();

        if (obr_.time().writeTime())
        {
            coordinateFrameState::write();
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::motionCoordinateFrame::frameVelocity
(
    const volVectorField& positions,
    bool addParentFrames
) const
{
    const vector omega = Omega();
    const vector vel = velocity().first();

    // Relative frame velocity
    tmp<volVectorField> Urf
    (
        new volVectorField
        (
            IOobject
            (
                "Urf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (
                (
                    dimensionedVector("Omega", dimless/dimTime, omega)
                 ^ (
                        positions
                      - dimensionedVector("CofR", dimLength, CofR())
                    )
                )
              + dimensionedVector("Ulinear", dimVelocity, vel)
            )
        )
    );

    volVectorField& Urfref(Urf.ref());

    if (validParentFrame() && addParentFrames)
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(parentFrameName_);
        Urfref += parentFrame.frameVelocity(positions, addParentFrames);
    }

    return Urf;
}


Foam::tmp<Foam::surfaceVectorField> Foam::motionCoordinateFrame::frameVelocity
(
    const surfaceVectorField& positions,
    bool addParentFrames
) const
{
    const vector omega = Omega();
    const vector vel = velocity().first();

    tmp<surfaceVectorField> Urf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "phirf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (
                (
                    dimensionedVector("Omega", dimless/dimTime, omega)
                  ^ (
                        positions
                      - dimensionedVector("CofR", dimLength, CofR())
                    )
                )
                + dimensionedVector("Ulinear", dimVelocity, vel)
            )
        )
    );
    surfaceVectorField& Urfref(Urf.ref());

    if (validParentFrame() && addParentFrames)
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(parentFrameName_);
        Urfref += parentFrame.frameVelocity(positions, addParentFrames);
    }

    return Urf;
}


Foam::tmp<Foam::vectorField> Foam::motionCoordinateFrame::frameVelocity
(
    const vectorField& positions,
    bool addParentFrames
) const
{
    tmp<vectorField> Urf =
        (Omega() ^ (positions - CofR())) + velocity().first();

    vectorField& Urfref = Urf.ref();

    if (validParentFrame() && addParentFrames)
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(parentFrameName_);
        Urfref += parentFrame.frameVelocity(positions, addParentFrames);
    }

    return Urf;
}


Foam::vector Foam::motionCoordinateFrame::frameVelocity
(
    const vector& position,
    bool addParentFrames
) const
{
    vector Urf = (Omega() ^ (position - CofR())) + velocity().first();

    if (validParentFrame() && addParentFrames)
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(parentFrameName_);
        Urf += parentFrame.frameVelocity(position, addParentFrames);
    }

    return Urf;
}


bool Foam::motionCoordinateFrame::writeData(Ostream& os) const
{
    os.writeEntry("origin", coorSys0().origin());
    os.writeEntry("axis", coorSys0().e3());
    return true;
}


// ************************************************************************* //
