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
    (c) 2011 OpenFOAM Foundation
    (c) 2013-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "pointPatchFields/derived/bufferZoneFixedDisplacement/bufferZoneFixedDisplacementPointPatchVectorField.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "patchDist/patchDistWave/patchDistWave.H"
#include "patchDist/wallPoint/wallPoint.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    radialVelocity_(vector::zero),
    linearVelocity_(vector::zero),
    CofG_(vector::zero),
    bufferSize_(0.0),
    p0_(p.localPoints())
{}


bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    radialVelocity_(dict.lookup("radialVelocity")),
    linearVelocity_(dict.lookup("linearVelocity")),
    CofG_(dict.lookup("CofG")),
    bufferSize_(readScalar(dict.lookup("bufferSize")))
{
    calculateBufferCellsAndPoints();

    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const bufferZoneFixedDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    radialVelocity_(ptf.radialVelocity_),
    linearVelocity_(ptf.linearVelocity_),
    CofG_(ptf.CofG_),
    bufferSize_(ptf.bufferSize_),
    p0_(mapper(ptf.p0_))
{}


bufferZoneFixedDisplacementPointPatchVectorField::
bufferZoneFixedDisplacementPointPatchVectorField
(
    const bufferZoneFixedDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    radialVelocity_(ptf.radialVelocity_),
    linearVelocity_(ptf.linearVelocity_),
    CofG_(ptf.CofG_),
    bufferSize_(ptf.bufferSize_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void bufferZoneFixedDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    m(p0_, p0_);
}


void bufferZoneFixedDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const bufferZoneFixedDisplacementPointPatchVectorField& aODptf =
        refCast<const bufferZoneFixedDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void bufferZoneFixedDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    vector transVec = linearVelocity_*t.value();

    vector eulerAngles = radialVelocity_*t.value();
    eulerAngles *= Foam::constant::mathematical::pi/180.0;

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-CofG_ + -transVec)*R*septernion(CofG_));

    vectorField lp = this->patch().localPoints();

    vectorField::operator=
    (
        transformPoints
        (
            TR,
            p0_
         )  - p0_
    );
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void bufferZoneFixedDisplacementPointPatchVectorField::manipulateMatrix
(
    const Field<vector>& oldPoints,
    labelList& bufferCells,
    Field<vector>& displacement
)
{
    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    vector transVec = linearVelocity_*t.value();

    vector eulerAngles = radialVelocity_*t.value();
    eulerAngles *= Foam::constant::mathematical::pi/180.0;

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-CofG_ + -transVec)*R*septernion(CofG_));

    bufferCells = bufferCells_;

    vectorField oldCellCentres(bufferCells.size(), vector::zero);
    forAll(bufferCells, i)
    {
        label cellI = bufferCells[i];
        cell c = mesh.cells()[cellI];
        oldCellCentres[i] = c.centre(oldPoints, mesh.faces());
    }

    displacement = transformPoints
    (
        TR,
        oldCellCentres
     )  -  oldCellCentres;
}


void bufferZoneFixedDisplacementPointPatchVectorField::setField
(
    const Field<vector>& oldPoints,
    Field<vector>& displacement
)
{
    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    vector transVec = linearVelocity_*t.value();

    vector eulerAngles = radialVelocity_*t.value();
    eulerAngles *= Foam::constant::mathematical::pi/180.0;

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-CofG_ + -transVec)*R*septernion(CofG_));

    UIndirectList<point>(displacement, bufferPoints_) = transformPoints
    (
        TR,
        pointField(oldPoints, bufferPoints_)
     )  -  pointField(oldPoints, bufferPoints_);
}


void Foam::bufferZoneFixedDisplacementPointPatchVectorField::calculateBufferCellsAndPoints()
{
    const polyMesh& mesh = this->internalField().mesh()();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    label patchI = bMesh.findPatchID(this->patch().name());

    if (patchI)
    {
        labelHashSet patchSet(1);
        patchSet.insert(patchI);

        scalarField y(mesh.nCells());
        patchDistWave::wave<wallPoint>(mesh, patchSet, y, false);

        boolList isBufferPt(mesh.nPoints(), true);

        DynamicList<label> cells(mesh.nCells());
        DynamicList<label> points(mesh.nCells());
        forAll(y, cellI)
        {
            if (y[cellI] < bufferSize_)
            {
                cells.append(cellI);
            }
            else
            {
                cell c = mesh.cells()[cellI];
                labelList meshPoints = c.labels(mesh.faces());
                forAll(meshPoints, ptI)
                {
                    isBufferPt[meshPoints[ptI]] = false;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            isBufferPt,
            orEqOp<bool>(),
            false
         );

        forAll(isBufferPt, ptI)
        {
            if (isBufferPt[ptI])
            {
               points.append(ptI);
            }
        }
        cells.shrink();
        points.shrink();

        bufferCells_ = cells;
        bufferPoints_ = points;
    }
    else
    {
        bufferCells_ = labelList();
        bufferPoints_ = labelList();
    }
}


void bufferZoneFixedDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("radialVelocity", radialVelocity_);
    os.writeEntry("CofG", CofG_);
    os.writeEntry("bufferSize", bufferSize_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    bufferZoneFixedDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
