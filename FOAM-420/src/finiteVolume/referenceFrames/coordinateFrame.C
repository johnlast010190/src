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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/coordinateFrame.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(coordinateFrame, 0);
    defineRunTimeSelectionTable(coordinateFrame, dictionary);
    addToRunTimeSelectionTable(coordinateFrame, coordinateFrame, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateFrame::coordinateFrame
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& frameName
)
:
    MeshObject<fvMesh, UpdateableMeshObject, coordinateFrame>(mesh, frameName),
    coordinateFrameState(mesh.thisDb(), dict, frameName),
    coorFrameReg_(mesh),
    frameDict_(dict),
    parentFrameName_
    (
        frameDict_.lookupOrAddDefault<word>("parentFrameName", word::null)
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::coordinateFrame::~coordinateFrame()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordinateFrame::initParents() const
{
    if (parents_.empty())
    {
        UPtrList<coordinateFrame> parents;
        label nNestedFrames = 0;
        bool nextFrame = true;
        const coordinateFrame* framePtr = this;
        while (nextFrame)
        {
            ++nNestedFrames;
            parents.setSize(nNestedFrames);
            parents.set(nNestedFrames - 1, const_cast<coordinateFrame*>(framePtr));
            if (framePtr->validParentFrame())
            {
                framePtr = &framePtr->parentFrame();
            }
            else
            {
                nextFrame = false;
            }
        }

        parents_.setSize(nNestedFrames);
        label n = 0;
        for (label i = parents.size() - 1; i >= 0; --i)
        {
            parents_.set(n, &parents[i]);
            ++n;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateFrame::localCoordSystem() const
{
    if
    (
        validParentFrame()
     && coorSys().coordinateTransforms().isLocalOnRead()
    && !coorSys().coordinateTransforms().isLocal()
    )
    {
        coordinateSystem localCS(parents().last().coorSys().clone());
        const label parentIndex = parents().size() - 2;

        if (parents()[parentIndex].coorSys().coordinateTransforms().isLocal())
        {
            FatalErrorInFunction
                << "Assumes that all frames are in global coordinates"
                << nl << abort(FatalError);
        }

        frameData::transformToLocal(localCS, parents()[parentIndex].coorSys());

        return localCS.clone();
    }
    else
    {
        return coorSys().clone();
    }
}


bool Foam::coordinateFrame::validParentFrame() const
{
    if (mesh_.foundObject<coordinateFrame>(parentFrameName_))
    {
        parentFramePtr_ =
            mesh_.lookupObjectRefPtr<coordinateFrame>(parentFrameName_);
        return true;
    }
    else if (parentFrameName_ != word::null)
    {
        parentFramePtr_ = &coordinateFrame::New(mesh_, parentFrameName_);
        return true;
    }
    return false;
}


const Foam::septernion& Foam::coordinateFrame::transformation
(
    bool addParentFrames,
    label nCorr
) const
{
    globalTransform_ = septernion::I;

    if (addParentFrames)
    {
        forAll(parents(), i)
        {
            if (parents()[i].isDynamic())
            {
                globalTransform_ *=
                    parents()[i].decoupledTransformation(nCorr);
            }
        }
    }
    else
    {
        if (isDynamic())
        {
            globalTransform_ = decoupledTransformation(nCorr);
        }
    }

    return globalTransform_;
}


void Foam::coordinateFrame::updateState() const
{
    if (!isUpdated() && validParentFrame())
    {
        // Request construction of old times
        oldTime().oldTime();

        // First store old time data, then update
        storeOldTimes();

        // Frame can get updated based on parent frame
        updateCoordinateSystem();

        // Update transformation, velocities and accelerations
        const label nCorr = nFrameCorrector() + 1;
        transformations_ = List<septernion>(nCorr, septernion::I);
        velocities_ = List<vectorTuple>(nCorr, vectorTuple(Zero, Zero));
        accelerations_ = List<vectorTuple>(nCorr, vectorTuple(Zero, Zero));

        updateIndex_ = obr_.time().timeIndex();
        if (obr_.time().writeTime())
        {
            coordinateFrameState::write();
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::coordinateFrame::frameVelocity
(
    const volVectorField& positions,
    bool addParentFrames
) const
{
    updateState();

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
            mesh_,
            dimensionedVector("vector", dimVelocity, Zero)
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


Foam::tmp<Foam::surfaceVectorField> Foam::coordinateFrame::frameVelocity
(
    const surfaceVectorField& positions,
    bool addParentFrames
) const
{
    updateState();

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
                    dimensionedVector("Omega", dimless/dimTime, Zero)
                  ^ (
                        positions
                      - dimensionedVector("orig", dimLength, Zero)
                    )
                )
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


Foam::tmp<Foam::vectorField> Foam::coordinateFrame::frameVelocity
(
    const vectorField& positions,
    bool addParentFrames
) const
{
    updateState();

    tmp<vectorField> Urf(positions*0);
    vectorField& Urfref(Urf.ref());

    if (validParentFrame() && addParentFrames)
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(parentFrameName_);
        Urfref += parentFrame.frameVelocity(positions, addParentFrames);
    }

    return Urf;
}


Foam::vector Foam::coordinateFrame::frameVelocity
(
    const vector& position,
    bool addParentFrames
) const
{
    updateState();

    if (validParentFrame() && addParentFrames)
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(parentFrameName_);
        return parentFrame.frameVelocity(position, addParentFrames);
    }

    return Zero;
}


void Foam::coordinateFrame::attachPatch(const label& patchi) const
{
    coorFrameReg_.attachPatch(patchi);

    if (validParentFrame())
    {
        const coordinateFrame& parentFrame =
            mesh_.lookupObject<coordinateFrame>(parentFrameName_);
        parentFrame.attachPatch(patchi);
    }
}


void Foam::coordinateFrame::inletFluxVelocity
(
    vectorField& Up,
    const label patchi,
    const objectRegistry& obr,
    const word& UName,
    bool addedMeshPhi
) const
{
    const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];

    // For nested dynamic frames remove normal component of velocity for all
    // frames but add meshPhi normal velocity component only ones
    if (mesh_.changing() && isDynamic())
    {
        const vectorField& nf = mesh_.boundaryMesh()[patchi].faceNormals();
        tmp<vectorField> Uframe = frameVelocity(Cf, false);
        Uframe.ref() -= nf*(nf & Uframe());

        // Mesh moving - normal component of velocity is removed and replaced
        // with the mesh flux velocity
        if (!addedMeshPhi)
        {
            const volVectorField& U = obr.lookupObject<volVectorField>(UName);
            const fvPatch& pp = U.boundaryField()[patchi].patch();
            scalarField phip
            (
                pp.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
            );
            tmp<scalarField> Un = phip/(pp.magSf() + VSMALL);
            Uframe.ref() += nf*Un;
            Up += Uframe;
            addedMeshPhi = true;
        }
    }
    else
    {
        Up += frameVelocity(Cf, false);
    }

    if (validParentFrame())
    {
        parentFramePtr_->inletFluxVelocity
        (
            Up,
            patchi,
            obr,
            UName,
            addedMeshPhi
        );
    }
}


void Foam::coordinateFrame::calculateBoundaryVelocity
(
    vectorField& Up,
    const label patchi,
    const objectRegistry& obr,
    const word& UName,
    const bool inletFlux
) const
{
    const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];
    if (inletFlux)
    {
        inletFluxVelocity(Up, patchi, obr, UName);
        return;
    }
    else if (!isAttachToMRF(patchi))
    {
        attachPatch(patchi);
    }
    const polyPatch& pp = mesh_.boundaryMesh()[patchi];
    const vectorField& nf = mesh_.boundaryMesh()[patchi].faceNormals();
    vectorField Uframe(frameVelocity(Cf, false));
    boolList isFaceInMRF(Up.size(), false);

    forAll(coorFrameReg_.registeredNames(), objI)
    {
        const word& objName = coorFrameReg_.registeredNames()[objI];
        if (mesh_.foundObject<frameSourceFaces>(objName))
        {
            const frameSourceFaces& frameSrcFaces =
                mesh_.lookupObject<frameSourceFaces>(objName);
            forAll(frameSrcFaces.includedFaces()[patchi], i)
            {
                const label patchFacei =
                    frameSrcFaces.includedFaces()[patchi][i];
                Up[patchFacei] += Uframe[patchFacei];
                isFaceInMRF[patchFacei] = true;
            }
        }
    }

    forAll(pp, patchFacei)
    {
        if (!isFaceInMRF[patchFacei])
        {
            const vector& faceNormal = nf[patchFacei];
            Up[patchFacei] +=
                Uframe[patchFacei]
              - faceNormal*(faceNormal & Uframe[patchFacei]);
        }
    }

    // If the parent frame isn't valid then it is the last frame in the chain
    // and the mesh flux velocity is added
    if (validParentFrame())
    {
        parentFramePtr_->calculateBoundaryVelocity
        (
            Up,
            patchi,
            obr,
            UName,
            inletFlux
        );
    }
    else if (mesh_.changing())
    {
        // Mesh moving - normal component of velocity is removed and replaced
        // with the mesh flux velocity
        const volVectorField& U = obr.lookupObject<volVectorField>(UName);
        const fvPatch& pp = U.boundaryField()[patchi].patch();
        scalarField phip
        (
            pp.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
        );
        tmp<scalarField> Un = phip/(pp.magSf() + VSMALL);
        Up += nf*Un;
    }
}


bool Foam::coordinateFrame::isAttachToMRF(const label& patchi) const
{
    if (coorFrameReg_.isAttachToMRF(patchi))
    {
        return true;
    }
    if (validParentFrame())
    {
        return parentFrame().isAttachToMRF(patchi);
    }
    return false;
}


// ************************************************************************* //
