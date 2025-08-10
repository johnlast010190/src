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
    (c) 2015-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/deformingBodyGIBFvMesh/deformingBodyGIBFvMesh.H"
#include "dynamicGIBFvMesh/deformingBodyMotionFunction/deformingBodyMotionFunction.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/volFields/volFields.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"
#include "triSurface/orientedSurface/orientedSurface.H"
#include "regionSplit/regionSplit.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "meshes/pointMesh/pointMesh.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "cfdTools/general/include/fvCFD.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "dynamicGIBFvMesh/movingGIBTools/twoDPointGIBCorrector/twoDPointGIBCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deformingBodyGIBFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        deformingBodyGIBFvMesh,
        IOobject
    );

    addToRunTimeSelectionTable
    (
        dynamicGIBFvMesh,
        deformingBodyGIBFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::deformingBodyGIBFvMesh::initialization()
{
    writeInterface_ =
        dynamicMeshCoeffs_.lookupOrDefault<Switch>("writeInterface", false);

    if (dynamicMeshCoeffs_.found("motionFunctions"))
    {
        const dictionary& mfDict =
            dynamicMeshCoeffs_.subDict("motionFunctions");
        DBMFPtr_ =
            deformingBodyMotionFunction::New(*this, mfDict, this->time());
    }

    this->faceZones().instance() = time().timeName();
    this->cellZones().instance() = time().timeName();

    postPro();

    Info<< dynamicMeshCoeffs_ << endl;
}


void Foam::deformingBodyGIBFvMesh::makeRecAllPoints0() const
{
    if (recAllPoints0Ptr_)
    {
        FatalErrorInFunction << abort(FatalError);
    }

    recAllPoints0Ptr_ = new pointField(*basePoints_);

    pointField& recAllPoints0 = *recAllPoints0Ptr_;
    tmp<vectorField> tOldBoundary = oldBoundaryLocation();
    recAllPoints0 = tOldBoundary.ref();
    syncPoints(recAllPoints0);
}


tmp<vectorField> Foam::deformingBodyGIBFvMesh::movePolyPatch
(
    primitivePatch& pp,
    const vectorField& interSpeed
)
{
    PrimitivePatchInterpolation<primitivePatch> pInterC(pp);

    scalarField interSpeedmag(mag(interSpeed));
    tmp<vectorField> pinterSpeedtmp
    (
        pInterC.faceToPointInterpolate(interSpeed)
    );
    tmp<scalarField> interSpeedmagtmp
    (
        pInterC.faceToPointInterpolate(interSpeedmag)
    );
    vectorField& pinterSpeed = pinterSpeedtmp.ref();
    pinterSpeed /= (mag(pinterSpeed) + SMALL);
    scalarField& pinterSpeedmag = interSpeedmagtmp.ref();

    vectorField pnf(pinterSpeed*pinterSpeedmag);

    labelList dispCount(pp.nPoints(), 1);

    syncTools::syncPointList
    (
        *this,
        pp.meshPoints(),
        pnf,
        plusEqOp<point>(),
        vector::zero,
        mapDistribute::transform()
    );
    syncTools::syncPointList
    (
        *this,
        pp.meshPoints(),
        dispCount,
        plusEqOp<label>(),
        label(0),
        mapDistribute::transform()
    );

    forAll(pnf, pI)
    {
        pnf[pI] = pnf[pI]/dispCount[pI];
    }

    tmp<vectorField> newPpPointst(new vectorField(this->points()));
    pointField& newPpPoints = newPpPointst.ref();

    const labelList& patchPoints = pp.meshPoints();
    forAll(patchPoints, pI)
    {
        const label& ppI = patchPoints[pI];
        newPpPoints[ppI] += time().deltaTValue()*pnf[pI];
    }
    return newPpPointst;
}


void Foam::deformingBodyGIBFvMesh::writeInterface
(
    const indirectPolyPatch& gibPolyPatch
)
{
    if (writeInterface_ && time().outputTime())
    {
        simpleVTKWriter
        (
            gibPolyPatch.localFaces(),
            gibPolyPatch.localPoints()
        ).write
        (
            postProFolder_
          + "/"
          + word("vis")
          + this->time().timeName()
          + ".vtk"
        );
    }
}


void Foam::deformingBodyGIBFvMesh::postPro()
{
    if (Pstream::master() || !Pstream::parRun())
    {
        postProFolder_ =
        (
            Pstream::parRun()
          ? this->time().path()/".."
          : this->time().path()
        );
        postProFolder_ += ("/postProcessing/GIBInterface");
        mkDir(postProFolder_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deformingBodyGIBFvMesh::deformingBodyGIBFvMesh
(
    const IOobject& io,
    const word& typeN
)
:
    dynamicGIBFvMesh(io, typeN),
    postProFolder_(),
    writeInterface_(false),
    DBMFPtr_(nullptr)
{
    deformingBodyGIBFvMesh::initialization();
}


Foam::deformingBodyGIBFvMesh::deformingBodyGIBFvMesh
(
    const IOobject& io
)
:
    dynamicGIBFvMesh(io, typeName),
    postProFolder_(),
    writeInterface_(false),
    DBMFPtr_(nullptr)
{
    deformingBodyGIBFvMesh::initialization();
}


Foam::deformingBodyGIBFvMesh::deformingBodyGIBFvMesh
(
    const IOobject& io,
    const dictionary dict
)
:
    dynamicGIBFvMesh(io, dict),
    postProFolder_(),
    writeInterface_(false),
    DBMFPtr_(nullptr)
{
    deformingBodyGIBFvMesh::initialization();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::deformingBodyGIBFvMesh::~deformingBodyGIBFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::deformingBodyGIBFvMesh::updateInit(const word& zoneName)
{
    IOobject triIO
    (
        triName_,
        this->time().constant(),
        "triSurface",
        *this,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    const fileName gibSurf = this->time().constant()/"triSurface"/triName_;
    if (isFile(gibSurf))
    {
        ibMeshPtr_ = new triSurfaceMesh(triIO);
    }
    else
    {
        FatalErrorInFunction
            << "Could not find file " << gibSurf << nl << exit(FatalError);
    }

    dynamicGIBFvMesh::updateInit(zoneName);
}


bool Foam::deformingBodyGIBFvMesh::update()
{
    storeOldTimes();
    oldPoints_ = points();
    prevPoints_ = points();
    DBMFPtr_->update();

    const fvPatch& gibPatch(this->boundary()[masterGIB_]);

    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(gibPatch.patch());

    faceList faces(gibPolyPatch);
    const boolList& flipMap = gibPolyPatch.fm();

    forAll(faces, fI)
    {
        if (flipMap[fI])
        {
            faces[fI].flip();
        }
    }
    pointField pointsF(this->points());
    primitivePatch pp(SubList<face>(faces, faces.size()), pointsF);

    pointsF = movePolyPatch(pp, DBMFPtr_->interfaceVelocity(masterGIB_)());

    //- important update pp (cleaning geometry)
    pp.clearGeom();

    mapGIB mapCl(*this, pp);

    //- !careful if this is destroyed ibMeshPtr_ reference is wrong
    deleteDemandDrivenData(ibMeshPtr_);
    ibMeshPtr_ =
        new triSurfaceMesh
        (
            IOobject
            (
                triName_,
                this->time().constant(),
                "triSurface",
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mapCl.triS()
        );

    writeInterface(gibPolyPatch);
    clearOutGIBData();
    doUpdate(mapCl, false);

    return true;
}


tmp<Foam::vectorField> Foam::deformingBodyGIBFvMesh::velocityCorrect
(
    const vectorField& pc
) const
{
    return DBMFPtr_->boundaryVelocity(masterGIB_);
}


tmp<Foam::vectorField>
Foam::deformingBodyGIBFvMesh::oldBoundaryLocation() const
{
    return DBMFPtr_->interfacePointsVelocity(masterGIB_);
}

// ************************************************************************* //
