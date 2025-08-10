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
    (c) 2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sweptCellMRF.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::sweptCellMRF::calculateConvectionField()
{
    boolList isAnMRFPatch(mesh_.boundary().size(), false);
    forAll(MRFPatchNames_, pI)
    {
        MRFPatchLabels_[pI] = mesh_.boundary().findPatchID(MRFPatchNames_[pI]);
        const label& pL = MRFPatchLabels_[pI];
        isAnMRFPatch[pL] = true;
        if (pL==-1)
        {
            FatalErrorInFunction
                << "Patch "<< MRFPatchNames_[pI] << " does not exist"
                << exit(FatalError);
        }
        fvPatchField<scalar>& boundaryConv = tracerField_.boundaryFieldRef()[pL];
        inletOutletFvPatchField<scalar>& mixedConv =
                dynamic_cast<inletOutletFvPatchField<scalar>&>(boundaryConv);
        mixedConv.refValue() = 1.;
        mixedConv.phiName() = "phiMRF";
    }

    surfaceScalarField phiMRF
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
    );

    phiMRF =
    (
        (axisOfRotation_ ^ (mesh_.Cf() - dimensionedVector("origin", dimLength, CofR_))) & mesh_.Sf()
    );

    surfaceScalarField::Boundary& phiMRFbf = phiMRF.boundaryFieldRef();
    forAll(phiMRFbf, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        if (!pp.coupled())
        {
            forAll(phiMRFbf[patchi], facei)
            {
                vector r = mesh_.Cf().boundaryField()[patchi][facei] - CofR_;
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

    setFvSchemesDict();
    setFvSolutionDict();

    for (int i=0; i<numberOfIterations_; i++)
    {
        fvScalarMatrix convEqn
        (
            fvm::div(phiMRF, tracerField_, "div(phiMRF,tracerField)")
        );
        forAll(mesh_.cells(), ci)
        {
            scalarField& diag = convEqn.diag();

            if (diag[ci]==0)
            {
                diag[ci]=1.;
            }
        }
        convEqn.relax();
        convEqn.solve();
    }

    volScalarField tracerField_1 = tracerField_;
    phiMRF = -phiMRF;
    tracerField_ *= 0;
    for (int i=0; i<numberOfIterations_; i++)
    {
        fvScalarMatrix convEqn
        (
            fvm::div(phiMRF, tracerField_, "div(phiMRF,tracerField)")
        );
        forAll(mesh_.cells(), ci)
        {
            scalarField& diag = convEqn.diag();

            if (diag[ci]==0)
            {
                diag[ci]=1.;
            }
        }
        convEqn.relax();
        convEqn.solve();
    }

    tracerField_ = max(tracerField_1, tracerField_);
}

void Foam::sweptCellMRF::setFvSchemesDict()
{
    fvMesh& refMesh = const_cast<fvMesh&>(mesh_);
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

    refMesh.schemes().setLocalSchemeDict(fvSchemesDict);
}

void Foam::sweptCellMRF::setFvSolutionDict()
{
    fvMesh& refMesh = const_cast<fvMesh&>(mesh_);
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
    fvSolutionDict.subDict("relaxationFactors").subDict("equations").add
    (
        "tracerField", 0.99
    );
    refMesh.solution().setLocalSolutionDict(fvSolutionDict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sweptCellMRF::sweptCellMRF
(
    const dictionary& MRFDict,
    const fvMesh& mesh
)
:
    MRFdict_(MRFDict),
    mesh_(mesh),
    MRFPatchNames_(MRFDict.lookup<wordList>("patches")),
    MRFPatchLabels_(MRFPatchNames_.size()),
    tracerField_
    (
        IOobject
        (
            "tracerField",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("scalar", dimless, 0),
        "inletOutlet"
    ),
    coorFramePtr_(coordinateFrame::lookupNew(mesh_, MRFDict)),
    axisOfRotation_(coorFramePtr_->axis0()),
    CofR_(coorFramePtr_->CofR0()),
    tracerThreshold_(MRFDict.lookupOrDefault<scalar>("tracerThreshold", 0.9)),
    writeTracerField_(MRFDict.lookupOrDefault<bool>("writeTracerField", false)),
    numberOfIterations_(MRFDict.lookupOrDefault<label>("numberOfIterations", 30))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sweptCellMRF::~sweptCellMRF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

labelList Foam::sweptCellMRF::getMRFCells()
{
    calculateConvectionField();

    boolList isAnMRFCell(mesh_.nCells(), true);

    scalar maxDist(-1);
    forAll(MRFPatchLabels_, pI)
    {
        const label& pL = MRFPatchLabels_[pI];
        const polyPatch& pp = mesh_.boundaryMesh()[pL];
        forAll(pp, fI)
        {
            vector r = mesh_.Cf().boundaryField()[pL][fI] - CofR_;
            vector rotN = axisOfRotation_/mag(axisOfRotation_);
            vector axisD = r - r*(r&rotN);
            maxDist = max(maxDist, mag(axisD));
        }
    }

    reduce(maxDist, maxOp<scalar>());
    forAll(mesh_.cells(), cI)
    {
        if (tracerField_[cI] >= tracerThreshold_ && isAnMRFCell[cI])
        {
            isAnMRFCell[cI] = true;
        }
        else
        {
            isAnMRFCell[cI] = false;
            tracerField_[cI] = 0;
        }
        vector r = mesh_.C()[cI] - CofR_;
        vector rotN = axisOfRotation_/mag(axisOfRotation_);
        vector axisD = r - r*(r&rotN);
        if (mag(axisD) > maxDist)
        {
            isAnMRFCell[cI] = false;
            tracerField_[cI] = 0;
        }
    }
    forAll(MRFPatchLabels_, pI)
    {
        const label& pL = MRFPatchLabels_[pI];
        const polyPatch& pp = mesh_.boundaryMesh()[pL];
        const labelList& fCells = pp.faceCells();
        forAll(tracerField_.boundaryField()[pL], faceI)
        {
            vector r = mesh_.Cf().boundaryField()[pL][faceI] - CofR_;
            vector fN = mesh_.Sf().boundaryField()[pL][faceI];
            vector fNMag = fN/mag(fN);
            vector rotDir = axisOfRotation_^r;
            if (mag(rotDir) > ROOTVSMALL)
            {
                vector rotDirMag = rotDir/mag(rotDir);

                if (mag(rotDirMag&fNMag) > 0.7)
                {
                    const label& cL = fCells[faceI];
                    tracerField_[cL] = 1;
                    isAnMRFCell[cL] = true;
                }
            }
        }
    }
    label nCells = 0;
    forAll(mesh_.cells(), cI)
    {
        if (isAnMRFCell[cI])
        {
            nCells++;
        }
    }
    labelList MRFcells(nCells);
    nCells = 0;
    forAll(mesh_.cells(), cI)
    {
        if (isAnMRFCell[cI])
        {
            MRFcells[nCells] = cI;
            nCells++;
        }
    }
    if (writeTracerField_)
    {
        if (MRFdict_.found("tracerName"))
        {
            word tracerName = MRFdict_.lookup("tracerName");
            volScalarField outputField
            (
                IOobject
                (
                    tracerName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                tracerField_
            );
            outputField.write();
        }
        else
        {
            WarningInFunction
                << "writeTracerField option was set to true."
                << nl
                << "   But no tracerName was defined. "<<nl
                << "   The default tracerName will be used for the outputField"
                << endl;
            tracerField_.write();
        }
    }

    return MRFcells;
}


// ************************************************************************* //
