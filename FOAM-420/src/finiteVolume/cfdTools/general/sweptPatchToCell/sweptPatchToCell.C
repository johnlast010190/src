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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/sweptPatchToCell/sweptPatchToCell.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sweptPatchToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, sweptPatchToCell, word);
}


Foam::topoSetSource::addToUsageTable Foam::sweptPatchToCell::usage_
(
    sweptPatchToCell::typeName,
    "\n    Usage:  \n\n"
    "    Select all cells within the swept volume of a patch\n\n"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sweptPatchToCell::calculateTracerField()
{
    setFluxNameInBoundaryPatches();
    setBoundaryValueInSweptPatches();
    constructRotationalFlux();

    setFvSchemesDict();
    setFvSolutionDict();

    solveConvectionEquation();

    volScalarField tracerField_1 = tracerFieldPtr_();
    reverseFluxAndInitializeTracer();

    solveConvectionEquation();

    tracerFieldPtr_() = max(tracerField_1, tracerFieldPtr_());

    findCellsInSet();
}

void Foam::sweptPatchToCell::setFluxNameInBoundaryPatches()
{
    forAll(mesh_.boundary(), pI)
    {
        fvPatchField<scalar>& boundaryConv =
            tracerFieldPtr_().boundaryFieldRef()[pI];
        try
        {
            inletOutletFvPatchField<scalar>& mixedConv
            = dynamic_cast<inletOutletFvPatchField<scalar>&>(boundaryConv);

            mixedConv.phiName() = "phiMRF";
        }
        catch(const std::bad_cast& e)
        {
            // Cast failed
        }
    }
}

void Foam::sweptPatchToCell::setBoundaryValueInSweptPatches()
{
    forAll(patchNames_, pI)
    {
        patchLabels_[pI] = mesh_.boundary().findPatchID(patchNames_[pI]);
        const label& pL = patchLabels_[pI];
        if (pL==-1)
        {
            FatalErrorInFunction
            << "Patch "<< patchNames_[pI] << " does not exist"
            << exit(FatalError);
        }
        fvPatchField<scalar>& boundaryConv =
            tracerFieldPtr_().boundaryFieldRef()[pL];
        inletOutletFvPatchField<scalar>& mixedConv =
            dynamic_cast<inletOutletFvPatchField<scalar>&>(boundaryConv);
        mixedConv.refValue() = 1.;
    }
}

void Foam::sweptPatchToCell::constructRotationalFlux()
{
    phiMRFPtr_() =
        (axisOfRotation_ ^ (mesh_.Cf() - dimensionedVector("origin",
        dimLength, origin_))) & mesh_.Sf();

    surfaceScalarField::Boundary& phiMRFbf = phiMRFPtr_().boundaryFieldRef();
    forAll(phiMRFbf, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        if (!pp.coupled())
        {
            forAll(phiMRFbf[patchi], facei)
            {
                vector r = mesh_.Cf().boundaryField()[patchi][facei] - origin_;
                vector fN = mesh_.Sf().boundaryField()[patchi][facei];
                vector fNMag = fN/mag(fN);
                vector rotDir = axisOfRotation_^r;
                if (mag(rotDir) > ROOTVSMALL)
                {
                    vector rotDirMag = rotDir/mag(rotDir);
                    if (mag(rotDirMag&fNMag) < 0.1)
                    {
                        phiMRFbf[patchi][facei] = 0;
                    }
                }
            }
        }
    }
}

void Foam::sweptPatchToCell::setFvSchemesDict()
{
    dictionary fvSchemesDict = mesh_.schemes().localSchemeDict();

    char hdg[] = ("cellLimited Gauss linear 1");
    fvSchemesDict.add("gradSchemes", dictionary(), false);
    fvSchemesDict.subDict("gradSchemes").add
    (
        string("grad(tracerField)"),hdg
    );

    char hdd[] = ("Gauss QUICK grad(tracerField)");
    fvSchemesDict.add("divSchemes", dictionary(), false);
    fvSchemesDict.subDict("divSchemes").add
    (
        string("div(phiMRF,tracerField)"),hdd
    );

    fvSchemes& fvs = const_cast<Foam::fvSchemes&>(mesh_.schemes());
    fvs.setLocalSchemeDict(fvSchemesDict);
}

void Foam::sweptPatchToCell::setFvSolutionDict()
{
    dictionary fvSolutionDict = mesh_.solution().localSolutionDict();
    dictionary hds;
    hds.add("solver", "smoothSolver");
    hds.add("smoother", "symGaussSeidel");
    hds.add("nSweeps", "2");
    hds.add("tolerance", "1e-14");
    hds.add("relTol", "0.1");

    fvSolutionDict.add("solvers", dictionary(), false);
    fvSolutionDict.subDict("solvers").add("tracerField", hds, false);
    fvSolutionDict.add("relaxationFactors", dictionary(), false);
    fvSolutionDict.subDict("relaxationFactors")
        .add("equations", dictionary(), false);
    fvSolutionDict.subDict("relaxationFactors").
        subDict("equations").add("tracerField", 0.99);
    fvSolution& fvs = const_cast<Foam::fvSolution&>(mesh_.solution());
    fvs.setLocalSolutionDict(fvSolutionDict);
}

void Foam::sweptPatchToCell::solveConvectionEquation()
{
    for (int i=0; i<numberOfSolvingIterations_; i++)
    {
        fvScalarMatrix convEqn
        (
            fvm::div(phiMRFPtr_(), tracerFieldPtr_(), "div(phiMRF,tracerField)")
        );
        correctZeroDiagonalElements(convEqn);
        convEqn.relax();
        convEqn.solve();
    }
}

void Foam::sweptPatchToCell::correctZeroDiagonalElements
(
    fvScalarMatrix& eqnMatrix
)
{
    forAll(mesh_.cells(), ci)
    {
        scalarField& diag = eqnMatrix.diag();

        if (diag[ci]==0)
        {
            diag[ci]=1.;
        }
    }
}

void Foam::sweptPatchToCell::reverseFluxAndInitializeTracer()
{
    phiMRFPtr_() = -phiMRFPtr_();
    tracerFieldPtr_() *= 0;
}

void Foam::sweptPatchToCell::findCellsInSet()
{
    excludeCellsWithValuesSmallerThanThreshold();
    excludeCellsFurtherThanMaximumDistance();
    forceAddBoundaryCellsOfSweptPatches();
}

void Foam::sweptPatchToCell::excludeCellsWithValuesSmallerThanThreshold()
{
    forAll(mesh_.cells(), cI)
    {
        if (tracerFieldPtr_()[cI]<tracerThresholdValue_)
        {
            cellBelongsToSet_[cI] = false;
        }
    }
}

void Foam::sweptPatchToCell::excludeCellsFurtherThanMaximumDistance()
{
    scalar maximumDistance = getMaximumDistanceFromTheAxis();
    forAll(mesh_.cells(), cI)
    {
        vector r = mesh_.C()[cI] - origin_;
        vector rotN = axisOfRotation_/mag(axisOfRotation_);
        vector axisD = r - r*(r&rotN);
        if (mag(axisD)>maximumDistance)
        {
            cellBelongsToSet_[cI] = false;
        }
    }
}

Foam::scalar Foam::sweptPatchToCell::getMaximumDistanceFromTheAxis()
{
    scalar maxDist(-1);
    forAll(patchLabels_, pI)
    {
        const label& pL = patchLabels_[pI];
        const polyPatch& pp = mesh_.boundaryMesh()[pL];
        forAll(pp, fI)
        {
            vector r = mesh_.Cf().boundaryField()[pL][fI] - origin_;
            vector rotN = axisOfRotation_/mag(axisOfRotation_);
            vector distanceFromAxis = r - r*(r&rotN);
            maxDist = max(maxDist, mag(distanceFromAxis));
        }
    }

    reduce(maxDist, maxOp<scalar>());
    return maxDist;
}

void Foam::sweptPatchToCell::forceAddBoundaryCellsOfSweptPatches()
{
    forAll(patchLabels_, pI)
    {
        const label& pL = patchLabels_[pI];
        const polyPatch& pp = mesh_.boundaryMesh()[pL];
        const labelList& fCells = pp.faceCells();
        forAll(tracerFieldPtr_().boundaryField()[pL], faceI)
        {
            vector r = mesh_.Cf().boundaryField()[pL][faceI] - origin_;
            vector fN = mesh_.Sf().boundaryField()[pL][faceI];
            vector fNMag = fN/mag(fN);
            vector rotDir = axisOfRotation_^r;
            if (mag(rotDir) > ROOTVSMALL)
            {
                vector rotDirMag = rotDir/mag(rotDir);

                if (mag(rotDirMag&fNMag) > 0.7)
                {
                    const label& cL = fCells[faceI];
                    tracerFieldPtr_()[cL] = 1;
                    cellBelongsToSet_[cL] = true;
                }
            }
        }
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::sweptPatchToCell::sweptPatchToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    mesh_
    (
        Foam::IOobject
        (
            mesh.name(),
            mesh.time().timeName(),
            mesh.time(),
            Foam::IOobject::NO_READ
         ),
         xferCopy(mesh.points()),
         xferCopy(mesh.faces()),
         xferCopy(mesh.faceOwner()),
         xferCopy(mesh.faceNeighbour())
    ),
    patchNames_(dict.lookup<wordList>("patches")),
    patchLabels_(patchNames_.size()),
    cellBelongsToSet_(mesh.nCells(), true),
    axisOfRotation_(dict.lookup("axis")),
    origin_(dict.lookup("origin")),
    tracerThresholdValue_
    (
        dict.lookupOrDefault<scalar>("tracerThreshold", 0.9)
    ),
    numberOfSolvingIterations_
    (
        dict.lookupOrDefault<label>("numberOfIterations", 30)
    )
{
    // Add the fvMesh boundary patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    List<polyPatch*> p(patches.size());

    forAll(p, patchi)
    {
        p[patchi] = patches[patchi].clone(mesh.boundaryMesh()).ptr();
    }

    mesh_.addFvPatches(p);

    tracerFieldPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "tracerField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("scalar", dimless, 0),
            "inletOutlet"
        )
    );

    phiMRFPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiMRF",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("scalar", dimVol, 0.0)
         )
    );

    calculateTracerField();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::sweptPatchToCell::~sweptPatchToCell() = default;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::sweptPatchToCell::combine(topoSet& set, const bool add) const
{
    forAll(mesh().cells(), celli)
    {
        if (cellBelongsToSet_[celli])
        {
            addOrDelete(set, celli, add);
        }
    }
}


void Foam::sweptPatchToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells within the swept volume of patches: "
            << patchNames_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells within the swept volume of patches: "
            << patchNames_ << endl;

        combine(set, false);
    }
}
