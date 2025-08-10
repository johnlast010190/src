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
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshRefinement/meshRefinement.H"
#include "patchDist/wallPointData/wallPointData.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "meshRefinement/fieldData/fieldData.H"
#include "sets/topoSets/cellSet.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "regionSplit/regionSplit.H"
#include "refinementSurfaces/refinementSurfaces.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::autoPtr<Foam::regionSplit> Foam::meshRefinement::splitWrappedMesh
(
    const boolList& cutoffCells
) const
{
    // Analyse regions.
    boolList blockedFace(mesh_.nFaces(), false);

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (cutoffCells[own[faceI]] && !cutoffCells[nei[faceI]])
        {
            blockedFace[faceI] = true;
        }
        else if (!cutoffCells[own[faceI]] && cutoffCells[nei[faceI]])
        {
            blockedFace[faceI] = true;
        }
    }

    polyBoundaryMesh& boundary
        = const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());

    forAll(boundary, patchI)
    {
        const polyPatch& pp = boundary[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start()+i;
                blockedFace[faceI] = cutoffCells[own[faceI]];
            }
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, xorEqOp());

        // Set region per cell based on walking
    return   autoPtr<regionSplit>
    (
        new regionSplit
        (
            mesh_,
            blockedFace
        )
    );
}


void Foam::meshRefinement::setFvDicts() const
{
    fvSchemes& schemes = const_cast<fvSchemes&>(mesh_.schemes());
    fvSolution& solution = const_cast<fvSolution&>(mesh_.solution());
    schemes.resetReadOpt(IOobject::READ_IF_PRESENT);
    solution.resetReadOpt(IOobject::READ_IF_PRESENT);

    //Clear out any starting fvSchemes, fvSolution
    dictionary& fvSchemesDict  = const_cast<dictionary&>(mesh_.schemes().dict());

    fvSchemesDict.clear();
    dictionary& fvSolutionDict  = const_cast<dictionary&>(mesh_.solution().dict());
    fvSolutionDict.clear();

    dictionary fvSchemes = mesh_.schemes().localSchemeDict();
    {
        if (!fvSchemes.found(word("divSchemes")))
        {
            fvSchemes.add(word("divSchemes"), dictionary(), false);
        }
        else
        {
            fvSchemes.subDict("divSchemes").clear();
        }

        char defaultDiv[]
            = "Gauss linear";

        fvSchemes.subDict("divSchemes").add
        (
            word("default"),
            defaultDiv,
            true
        );
        char phiUDiv[]
            = "Gauss upwind";

        fvSchemes.subDict("divSchemes").add
        (
            word("upwind"),
            phiUDiv,
            true
        );

        if (!fvSchemes.found(word("gradSchemes")))
        {
            fvSchemes.add(word("gradSchemes"), dictionary(), false);
        }
        else
        {
            fvSchemes.subDict("gradSchemes").clear();
        }

        char defaultGrad[]
            = "Gauss linear";

        fvSchemes.subDict("gradSchemes").add
        (
            word("default"),
            defaultGrad,
            true
        );

        if (!fvSchemes.found(word("interpolationSchemes")))
        {
            fvSchemes.add(word("interpolationSchemes"), dictionary(), false);
        }
        else
        {
            fvSchemes.subDict("interpolationSchemes").clear();
        }

        char defaultInterpolation[]
            = "linear";

        fvSchemes.subDict("interpolationSchemes").add
        (
            word("default"),
            defaultInterpolation,
            true
        );

        if (!fvSchemes.found(word("laplacianSchemes")))
        {
            fvSchemes.add(word("laplacianSchemes"), dictionary(), false);
        }
        else
        {
            fvSchemes.subDict("laplacianSchemes").clear();
        }

        char defaultLaplacian[]
            = "Gauss linear limited 0.333";

        fvSchemes.subDict("laplacianSchemes").add
        (
            word("default"),
            defaultLaplacian,
            true
        );

        if (!fvSchemes.found(word("snGradSchemes")))
        {
            fvSchemes.add(word("snGradSchemes"), dictionary(), false);
        }
        else
        {
            fvSchemes.subDict("snGradSchemes").clear();
        }

        char defaultSnGrad[]
            = "limited 0.333";

        fvSchemes.subDict("snGradSchemes").add
        (
            word("default"),
            defaultSnGrad,
            true
        );

        if (!fvSchemes.found(word("fluxRequired")))
        {
            fvSchemes.add(word("fluxRequired"), dictionary(), false);
        }
        else
        {
            fvSchemes.subDict("fluxRequired").clear();
        }

        char defaultFlux[]
            = "no";

        fvSchemes.subDict("fluxRequired").add
        (
            word("default"),
            defaultFlux,
            true
        );

        char pFlux[]
            = "";

        fvSchemes.subDict("fluxRequired").add
        (
            word("p"),
            pFlux,
            true
        );
    }
    mesh_.schemes().setLocalSchemeDict(fvSchemes);

    dictionary fvSolution = mesh_.solution().localSolutionDict();
    {
        if (!fvSolution.found(word("solvers")))
        {
            fvSolution.add(word("solvers"), dictionary(), false);
        }
        else
        {
            fvSolution.subDict("solvers").clear();
        }

        dictionary solver;
        solver.add("solver", "GAMG");
        solver.add("agglomerator", "faceAreaPair");
        solver.add("mergeLevels", "1");
        solver.add("nCellsInCoarsestLevel", "200");
        solver.add("tolerance", "1e-8");
        solver.add("relTol", 0.001);
        solver.add("smoother", "GaussSeidel");
        solver.add("nPreSweeps", "0");
        solver.add("nPostSweeps", "2");
        solver.add("nFinestSweeps", "2");

        fvSolution.subDict("solvers").add
        (
            "p",
            solver,
            true
        );
    }
    mesh_.solution().setLocalSolutionDict(fvSolution);
}


void Foam::meshRefinement::sourceSetup
(
    const scalar& volDist,
    const scalar& Vs,
    volScalarField& dVs
) const
{
    DynamicList<label> changedFaces(mesh_.nFaces()/100 + 100);
    DynamicList<wallPointData<label>> faceDist(changedFaces.size());

    forAll(mesh_.boundaryMesh(), patchI)
    {
         const polyPatch& patch = mesh_.boundaryMesh()[patchI];

         if
         (
            (!(patch.physicalType()=="outlet") && !(patch.type()=="outlet" ))
            && !isA<processorPolyPatch>(patch)
         )
         {
              forAll(patch.faceCentres(), patchFaceI)
              {
                   label meshFaceI = patch.start() + patchFaceI;
                   changedFaces.append(meshFaceI);

                   faceDist.append
                   (
                        wallPointData<label>
                        (
                             patch.faceCentres()[patchFaceI],
                             patchI, // passive label
                             0.0
                        )
                   );
              }
         }
    }
    changedFaces.shrink();
    faceDist.shrink();

    // Perform mesh wave to calculate distance to patch
    List<wallPointData<label>> faceInfo(mesh_.nFaces());
    List<wallPointData<label>> cellInfo(mesh_.nCells());

    FaceCellWave<wallPointData<label>> waveInfo
    (
        mesh_,
        changedFaces,
        faceDist,
        faceInfo,
        cellInfo,
        mesh_.globalData().nTotalCells()    // max iterations
    );

    scalar volDistSqr =  volDist * volDist;

    forAll(mesh_.cells(), cellI)
    {
        scalar dSqr = cellInfo[cellI].distSqr();

        if (dSqr < volDistSqr)
        {
            dVs[cellI] =  Vs;
        }
    }
}


Foam::labelList Foam::meshRefinement::solve
(
    const dictionary& dict,
    const pointField& locationsInMesh,
    const labelList& keepCells
) const
{
    Info<< nl << "Calculating Euler flow" << endl;

    //check availability of fvSchemes and fvSolution settings
    setFvDicts();

    int nNonOrthCorr = 0;
    if (mesh_.solution().dict().found(word("SIMPLE")))
    {
        dictionary simple = mesh_.solution().dict().subDict("SIMPLE");
        nNonOrthCorr =
            simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    }

    scalar c1 = dict.lookupOrDefault("C1", 1);
    scalar c2 = dict.lookupOrDefault("C2", 4);
    scalar Vs = dict.lookupOrDefault("Vs", 0.001);
    scalar volDistance = dict.lookupOrDefault("volDistance", 1.0);

    Switch volSources = dict.lookupOrDefault("volSources", true);
    Switch meshedInMM = dict.lookupOrDefault("meshInMM", false);

    label maxIter = dict.lookupOrDefault("maxIter", 50);

    dimensionedScalar C1("C1", dimless/dimTime, c1);
    dimensionedScalar C2("C2", dimless/dimLength, c2);

    // If meshing in mm need to rescale damping and volume sources
    if (meshedInMM)
    {
        C1 *= 1e-3;
        C2 *= 1e-3;
        Vs *= 1e-3;
    }

    const fvPatchList& patches = mesh_.boundary();

    polyBoundaryMesh& boundary
        = const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());
    fvBoundaryMesh& fvPatches
        = const_cast<fvBoundaryMesh&>(mesh_.boundary());

    HashTable<word> outletPatches;
    if (dict.found("outlets"))
    {
        wordList outPatches = wordList(dict.lookup("outlets"));
        forAll(outPatches, nI)
        {
            if (!outletPatches.found(outPatches[nI]))
            {
                outletPatches.insert(outPatches[nI], outPatches[nI]);
            }
        }
    }
    else
    {
        labelHashSet meshPatches(meshedPatches());
        forAll(patches, patchI)
        {
            word bType = patches[patchI].type();
            word name = patches[patchI].name();
            if (bType != "processor" && !meshPatches.found(patchI))
            {
                if (!outletPatches.found(name))
                {
                    outletPatches.insert(name, name);

                }
            }
        }
    }

    label sz = 0;
    forAll(patches, patchI)
    {
        if (outletPatches.found(patches[patchI].name()))
        {
            sz += patches[patchI].size();
        }
    }

    if (!outletPatches.size())
    {
        WarningInFunction
            << "Wrapping turned off as no outlet flow patches found "
            << "for wrapping solver. "
            << "Available patches: "<< boundary.names() << endl;
        return labelList(0);
    }
    else if (!returnReduce(sz, sumOp<label>()))
    {
        WarningInFunction
            << "Wrapping turned off as specified oulet flow patches "
            << "are of zero size." << endl;
        return labelList(0);
    }
    else
    {
        Info<<"Setting the following patches to flow outlet patches: "
            <<outletPatches.toc()<<endl;
    }


    dictionary ufieldDict;
    dictionary pfieldDict;

    ufieldDict.add("internalField", "uniform (0 0 0)");
    ufieldDict.add("boundaryConditions", dictionary(), false);

    pfieldDict.add("internalField", "uniform 0");
    pfieldDict.add("boundaryConditions", dictionary(), false);

    dictionary& uBCDict =
        ufieldDict.subDict("boundaryConditions");
    dictionary& pBCDict =
        pfieldDict.subDict("boundaryConditions");

    wordList updatedType(patches.size(), word());

    forAll(patches, patchI)
    {
        const word& name = patches[patchI].name();

        word bType = patches[patchI].type();
        uBCDict.add(name, dictionary(), false);
        pBCDict.add(name, dictionary(), false);
        dictionary& uDict = uBCDict.subDict(name);
        dictionary& pDict = pBCDict.subDict(name);

        bool reset = false;

        if (outletPatches.found(patches[patchI].name()))
        {
            updatedType[patchI] = patches[patchI].type();
            bType = "outlet";
            reset = true;

            uDict.add("type", "inletOutlet");
            uDict.add("inletValue", "uniform (0 0 0)");
            uDict.add("value", "uniform (0 0 0)");
            pDict.add("type", "fixedValue");
            pDict.add("value", "uniform 0");
        }
        else
        {
            if (bType == "processor")
            {
                pDict.add("type", "processor");
                pDict.add("value", "uniform 0");
                uDict.add("type", "processor");
                uDict.add("value", "uniform (0 0 0)");
            }
            else
            {

                pDict.add("type", "zeroGradient");
                if (volSources)
                {
                    uDict.add("type", "fixedValue");
                    uDict.add("value", "uniform (0 0 0)");
                }
                else
                {
                    updatedType[patchI] = patches[patchI].type();
                    bType = "inlet";
                    reset = true;
                    uDict.add("type", "surfaceNormalFixedValue");
                    uDict.add("refValue", "uniform -1");
                }
            }
        }
        if (reset)
        {
            boundary.set
            (
                patchI,
                polyPatch::New
                (
                    bType,
                    name,
                    patches[patchI].patch().size(),  // size
                    patches[patchI].patch().start(), // start
                    patchI,
                    boundary
                 )
             );

            fvPatches.set
            (
                patchI,
                fvPatch::New
                (
                    boundary[patchI],
                    // point to newly added polyPatch
                    mesh_.boundary()
                 )
             );
        }
    }

    autoPtr<GeometricField<vector, fvPatchField, volMesh>> UPtr
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "U",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensioned<vector>
            (
                "U",
                dimensionSet(0, 1, -1, 0, 0),
                pTraits<vector>::zero
             )
         )
     );

    UPtr->primitiveFieldRef()
        = Field<vector>("internalField", ufieldDict, mesh_.nCells());

    UPtr->boundaryFieldRef().reset
    (
        GeometricField<vector, fvPatchField, volMesh>::Boundary
        (
            mesh_.boundary(),
            UPtr->internalField(),
            boundarySetup<vector>
            (
                mesh_,
                UPtr->internalField(),
                ufieldDict
             )
         )
     );

    autoPtr<GeometricField<scalar, fvPatchField, volMesh>> pPtr
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                "p",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensioned<scalar>
            (
                "p",
                dimensionSet(0, 2, -2, 0, 0),
                pTraits<scalar>::zero
             )
         )
     );

    pPtr->primitiveFieldRef()
        = Field<scalar>("internalField", pfieldDict, mesh_.nCells());

    pPtr->boundaryFieldRef().reset
    (
        GeometricField<scalar, fvPatchField, volMesh>::Boundary
        (
            mesh_.boundary(),
            pPtr->internalField(),
            boundarySetup<scalar>
            (
                mesh_,
                pPtr->internalField(),
                pfieldDict
             )
         )
     );

    volScalarField dVs
    (
         IOobject
         (
              "dVs",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false
         ),
         mesh_,
         dimensionedScalar("one", dimless/dimTime, 0),
         zeroGradientFvPatchScalarField::typeName
    );

    volScalarField solid
    (
         IOobject
         (
              "solid",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE,
              false
         ),
         mesh_,
         dimensionedScalar("one", dimless, 1),
         zeroGradientFvPatchScalarField::typeName
    );

    GeometricField<vector, fvPatchField, volMesh>& U = UPtr();
    GeometricField<scalar, fvPatchField, volMesh>& p = pPtr();

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U) & mesh_.Sf()
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;

    if (mesh_.solution().dict().found(word("SIMPLE")))
    {
        setRefCell
        (
            p,
            mesh_.solution().dict().subDict("SIMPLE"),
            pRefCell,
            pRefValue
         );
    }

    if (volSources)
    {
        sourceSetup(volDistance, Vs, dVs);
    }

    pointField excludePoints = pointField();
    if (dict.found("excludePoints"))
    {
        excludePoints = pointField(dict.lookup("excludePoints"));
    }

    scalar defaultCutoff = 0.02;
    if (!volSources)
    {
        defaultCutoff = 1e5;
    }
    scalar cutoff = dict.lookupOrDefault("cutoff", defaultCutoff);
    scalar sigma = dict.lookupOrDefault("sigma", 3.0);
    label minRegionSize = dict.lookupOrDefault("minRegionSize", 5);
    bool writeFields = dict.lookupOrDefault("writeFields", false);
    bool invert = dict.lookupOrDefault("invert", false);

    if (!writeFields)
    {
        //Turn off solver debug
        Foam::lduMatrix::debug = 0;
    }

    DynamicList<label> cellsToRemove(mesh_.nCells());
    boolList cutoffCells(mesh_.nCells(), false);
    boolList removedCells(mesh_.nCells(), false);

    bool solve = true;
    label nIter = 1;
    label checkStage = 20;

    //debug
    if (writeFields)
    {
        Info<< "Writing pre-wrap mesh " << endl;
        const_cast<Time&>(mesh_.time())++;
        U.write();
        p.write();
        dVs.write();
        mesh_.write();
    }

    while (solve)
    {
        if (nIter % 10 == 0)
        {
            Info<< "Time = " << nIter << nl << endl;
        }

        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
            tmp<fvVectorMatrix> UEqn
            (
                fvm::div(phi, U, "upwind")
                + fvm::SuSp(-fvc::div(phi), U)
                + fvm::Sp((C1 + C2*mag(U)) * solid, U)
            );

            volScalarField rUA("rUA", 1.0/UEqn().A());
            p.boundaryFieldRef().updateCoeffs();

            U = UEqn().H() * rUA;
            UEqn.clear();
            phi = fvc::interpolate(U) & mesh_.Sf();
            adjustPhi(phi, U, p);

            // Non-orthogonal pressure corrector loop
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi) - dVs
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = fvc::reconstruct(phi);

            U.correctBoundaryConditions();
        }

        Info<<"Iter: "<< nIter <<" Max U (wrapper): "
            <<gMax(mag(U.primitiveField()))
            <<" Max p (wrapper): "
            <<gMax(p.primitiveField())<<endl;

        forAll(locationsInMesh, locationi)
        {
            label keepcelli = keepCells[locationi];
            if (keepcelli != -1)
            {
                Pout<<"Wrapper pressure : " << p.primitiveField()[keepcelli]
                    << " at locationInMesh : " << locationsInMesh[locationi]
                    <<endl;
            }
        }

        label nRemoved = 0;
        if (nIter % checkStage == 0)
        {
            //After first check increase check frequency
            checkStage = 10;
            if (excludePoints.size() == 0)
            {
                cutoffCells = false;
                forAll(mesh_.cells(), cellI)
                {
                    if (p[cellI] > cutoff)
                    {
                        cutoffCells[cellI] = true;
                    }
                }
                regionSplit cellRegion = splitWrappedMesh(cutoffCells);

                label nRegions = cellRegion.nRegions();
                scalarField maxZonePressure(nRegions, -GREAT);
                scalarField aveZonePressure(nRegions, 0.0);
                scalarField regionVolume(nRegions, 0.0);
                vectorField regionCentre(nRegions, vector::zero);
                forAll(cellRegion, cellI)
                {
                    label regionI = cellRegion[cellI];
                    maxZonePressure[regionI] =
                        max(maxZonePressure[regionI], p[cellI]);
                    aveZonePressure[regionI] +=
                        p[cellI]*mesh_.cellVolumes()[cellI];
                    scalar vol = mesh_.cellVolumes()[cellI];
                    regionVolume[regionI] += vol;
                    regionCentre[regionI] += vol*mesh_.cellCentres()[cellI];
                }
                Pstream::listCombineReduce(maxZonePressure, maxOp<scalar>());
                Pstream::listCombineReduce(aveZonePressure, plusOp<scalar>());
                Pstream::listCombineReduce(regionVolume, plusOp<scalar>());
                Pstream::listCombineReduce(regionCentre, plusOp<vector>());

                forAll(aveZonePressure, regionI)
                {
                    aveZonePressure[regionI] /= (regionVolume[regionI] + SMALL);
                }

                scalarField regionStdDev(nRegions, 0.0);
                forAll(cellRegion, cellI)
                {
                    label regionI = cellRegion[cellI];
                    regionStdDev[regionI] +=
                        sqr(p[cellI]-aveZonePressure[regionI])
                        *mesh_.cellVolumes()[cellI];
                }
                Pstream::listCombineReduce(regionStdDev, plusOp<scalar>());

                Info<<"Wrapper Statistics: "<<nl<<endl;
                forAll(regionStdDev, regionI)
                {
                    regionStdDev[regionI] = sqrt
                        (regionStdDev[regionI]/(regionVolume[regionI] + SMALL));

                    Info<<nl;
                    if (maxZonePressure[regionI] > cutoff)
                    {
                        Info<<"Active Region: "<<regionI<<" volume: "
                            <<regionVolume[regionI]<<endl;
                        Info<<"Region Centre" << regionCentre[regionI]
                            /(regionVolume[regionI] + SMALL)<<endl;
                        Info<<"Ave "<<aveZonePressure[regionI]<<"  Max: "
                            <<maxZonePressure[regionI]
                            <<" stdDev: "<<regionStdDev[regionI]
                            <<" Cutoff: "<<aveZonePressure[regionI]
                            - sigma * regionStdDev[regionI]
                            <<nl<<endl;
                    }
                    else
                    {
                        Info<<"Inactive Region: "<<regionI<<" volume: "
                            <<regionVolume[regionI]<<endl;
                        Info<<"Region Centre" << regionCentre[regionI]
                            /(regionVolume[regionI] + SMALL)<<endl;
                        Info<<"Ave "<<aveZonePressure[regionI]<<"  Max: "
                            <<maxZonePressure[regionI]
                            <<" stdDev: "<<regionStdDev[regionI]
                            <<nl<<endl;
                    }
                }

                cutoffCells = false;
                forAll(cellRegion, cellI)
                {
                    label regionI = cellRegion[cellI];
                    scalar zonePressure = maxZonePressure[regionI];
                    if
                    (
                        zonePressure > cutoff
                        && p[cellI] > (aveZonePressure[regionI]
                                       - sigma * regionStdDev[regionI])
                    )
                    {
                        cutoffCells[cellI] =  true;
                    }
                }
                //Remove islands
                regionSplit updatedCellRegion = splitWrappedMesh(cutoffCells);
                label nUpdatedRegions = updatedCellRegion.nRegions();
                labelList regionSize(nUpdatedRegions, 0);
                forAll(cellRegion, cellI)
                {
                    label regionI = updatedCellRegion[cellI];
                    regionSize[regionI]++;
                }
                Pstream::listCombineReduce(regionSize, plusOp<label>());

                forAll(cellRegion, cellI)
                {
                    label regionI = updatedCellRegion[cellI];
                    if
                    (
                        !removedCells[cellI] && cutoffCells[cellI]
                        && regionSize[regionI] > minRegionSize
                    )
                    {
                        removedCells[cellI] = true;
                        nRemoved++;
                    }
                }

                Info<<"Iter: "<< nIter <<" Max U (wrapper): "
                    <<gMax(mag(U.primitiveField()))
                    <<" Max p (wrapper): "
                    <<gMax(p.primitiveField())<<endl;

                if (volSources)
                {
                    forAll(cutoffCells, cellI)
                    {
                        if (cutoffCells[cellI])
                        {
                            dVs[cellI] = 0.0;
                            //solid[cellI] = GREAT;
                        }
                    }
                }
                else
                {
                    solve = false;
                }

                //debug
                if (writeFields)
                {
                    const_cast<Time&>(mesh_.time())++;
                    Info<< "Writing wrapping fields to time: "
                         <<  mesh_.time().timeName() <<endl;
                    U.write();
                    p.write();
                    dVs.write();

                    volScalarField wrapRegion
                    (
                        IOobject
                        (
                            "wrapRegions",
                            timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE,
                            false
                        ),
                        mesh_,
                        dimensionedScalar("zero", dimless, -1),
                        zeroGradientFvPatchScalarField::typeName
                    );

                    forAll(mesh_.cells(), cellI)
                    {
                        wrapRegion[cellI] = scalar(cellRegion[cellI]);
                    }
                    wrapRegion.write();
                }

                if (!returnReduce(nRemoved, sumOp<label>()))
                {
                    solve = false;
                }
                else
                {
                    Info<<nl;
                    Info<<"Removing "<<returnReduce(nRemoved, sumOp<label>())
                        <<" cells from mesh after iteration "<<nIter<<endl;
                    Info<<nl;
                }
            }
            else
            {
                solve = false;
            }
        }

        nIter++;
        if (nIter > maxIter)
        {
            solve = false;
        }
    }

    if (excludePoints.size() > 0)
    {
        surfaceScalarField pf( fvc::interpolate(p) );
        List<fieldData> allCellInfo(mesh_.nCells());
        List<fieldData> allFaceInfo(mesh_.nFaces());
        // Labels of seed faces
        DynamicList<label> seedFaces(mesh_.nFaces());
        //  data on seed faces
        DynamicList<fieldData> seedFacesInfo(mesh_.nFaces());

        // Dummy additional info for FaceCellWave
        int dummyTrackData = 0;
        forAll(excludePoints, i)
        {
            label cellI = mesh_.findCell
            (
                excludePoints[i],
                polyMesh::FACE_PLANES
            );

            if (cellI != -1)
            {
                const cell& cFaces = mesh_.cells()[cellI];
                forAll(cFaces, cFaceI)
                {
                    label faceI = cFaces[cFaceI];
                    if (!allFaceInfo[faceI].valid(dummyTrackData))
                    {
                        seedFaces.append(faceI);
                        seedFacesInfo.append
                        (
                            fieldData
                            (
                                pf[faceI],
                                pf[faceI],
                                true
                             )
                         );
                        allFaceInfo[faceI] =
                            seedFacesInfo[seedFacesInfo.size()-1];
                    }
                }
            }
        }

        forAll(mesh_.cells(), cellI)
        {
            if (!allCellInfo[cellI].valid(dummyTrackData))
            {
                allCellInfo[cellI] = fieldData
                                     (
                                         p[cellI],
                                         p[cellI],
                                         false
                                      );
            }
        }

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            if (!allFaceInfo[faceI].valid(dummyTrackData))
            {
                allFaceInfo[faceI] = fieldData
                                     (
                                         pf[faceI],
                                         pf[faceI],
                                         false
                                      );
            }
        }
        forAll(p.boundaryField(), patchI)
        {
            const fvPatchScalarField& bField = p.boundaryField()[patchI];
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
            label startFaceI = patch.start();

            forAll(bField, faceI)
            {
                label meshFaceI = startFaceI + faceI;
                if (!allFaceInfo[meshFaceI].valid(dummyTrackData))
                {
                    scalar boundFaceP = bField[faceI];
                    allFaceInfo[meshFaceI] = fieldData
                                             (
                                                 boundFaceP,
                                                 boundFaceP,
                                                 false
                                              );
                }
            }
        }

        // face-cell-face transport engine
        FaceCellWave<fieldData, int> propagationCalc
        (
            mesh_,
            allFaceInfo,
            allCellInfo,
            dummyTrackData
         );

        propagationCalc.setFaceInfo
        (
            seedFaces.shrink(),
            seedFacesInfo.shrink()
         );
        seedFaces.clear();
        seedFacesInfo.clear();
        propagationCalc.iterate(mesh_.globalData().nTotalFaces());

        forAll(allCellInfo, cellI)
        {
            if (!invert && allCellInfo[cellI].s())
            {
                cellsToRemove.append(cellI);
            }
            else if (invert && !allCellInfo[cellI].s())
            {
                cellsToRemove.append(cellI);
            }
        }
    }
    else
    {
        forAll(removedCells, cellI)
        {
            if (!invert && removedCells[cellI])
            {
                cellsToRemove.append(cellI);
            }
            else if (invert && !removedCells[cellI])
            {
                cellsToRemove.append(cellI);
            }
        }
    }
    cellsToRemove.shrink();

    // reset to original type
    forAll(patches, patchI)
    {
        if (updatedType[patchI].size())
        {
            const word& name = patches[patchI].name();
            word bType = updatedType[patchI];
            boundary.set
            (
                patchI,
                polyPatch::New
                (
                    bType,
                    name,
                    patches[patchI].patch().size(),  // size
                    patches[patchI].patch().start(), // start
                    patchI,
                    boundary
                 )
             );

            fvPatches.set
            (
                patchI,
                fvPatch::New
                (
                    boundary[patchI],
                    // point to newly added polyPatch
                    mesh_.boundary()
                 )
             );
        }
    }

    return labelList(cellsToRemove, true);
}


void Foam::meshRefinement::splitWrappedMesh
(
    const refinementParameters& refineParams,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const dictionary& dict,
    const dictionary& motionDict,
    const bool handleSnapProblems,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const writer<scalar>& leakPathFormatter
)
{
    const dictionary& wrapDict = dict.subDict("wrapper");

    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const pointField& locationsOutsideMesh = refineParams.locationsOutsideMesh();

    labelList keepCells(locationsInMesh.size(), -1);
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];

        keepCells[i] = findCell
        (
            insidePoint,
            mesh_,
            meshCutter_
        );
    }

    labelList cellsToRemove = solve(wrapDict,locationsInMesh,keepCells);

    //Check for removal of locationsInMesh cells
    forAll(keepCells,i)
    {
        const point& insidePoint = locationsInMesh[i];
        label keepCellI = keepCells[i];
        bool removeKeepCell = false;
        forAll(cellsToRemove, cI)
        {
            label celli = cellsToRemove[cI];
            if (celli == keepCellI)
            {
                removeKeepCell = true;
            }
        }

        if (returnReduce(removeKeepCell, orOp<bool>()))
        {
            FatalErrorInFunction
                << "Wrapper will remove cells at locationInMesh : "
                << insidePoint << nl
                << "Meshing needs to be rerun with"
                << " 1) wrapper deactivated"
                << " 2) change to cutoff method"
                << " 3) wrap keyword writeFields set to true." << nl
                << exit(FatalError);
        }
    }

    // Remove cells
    removeCells cellRemover(mesh_);

    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    label nExposed = returnReduce(exposedFaces.size(), sumOp<label>());

    label newPatchID = -1;

    if (nExposed > 0)
    {
        // Add patch
        dictionary patchInfo;
        patchInfo.set("type", wallPolyPatch::typeName);
        newPatchID = addPatch(mesh_, "wrapped", patchInfo);
    }

    doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        labelList(exposedFaces.size(), newPatchID),
        cellRemover
    );

    labelList ownPatch, nbrPatch;
    getBafflePatches
    (
        refineParams,
        globalToMasterPatch,
        ownPatch,
        nbrPatch
    );

    createBaffles(ownPatch, nbrPatch);

    if (handleSnapProblems)
    {
        Info<< nl
            << "Introducing baffles to block off problem cells" << nl
            << "----------------------------------------------" << nl
            << endl;

        const labelList surfacesToBaffle
        (
            surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones())
        );

        labelList facePatch
        (
            markFacesOnProblemCells
            (
                motionDict,
                refineParams,
                removeEdgeConnectedCells,
                perpendicularAngle,
                globalToMasterPatch,
                meshedPatches(),
                surfacesToBaffle
            )
            //markFacesOnProblemCellsGeometric(motionDict)
        );

        // Create baffles with same owner and neighbour for now.
        createBaffles(facePatch, facePatch);
    }

    splitMeshRegions
    (
        globalToMasterPatch,
        globalToSlavePatch,
        locationsInMesh,
        locationsOutsideMesh,
        leakPathFormatter,
        refineParams.fullLeakChecks()
    );
}
