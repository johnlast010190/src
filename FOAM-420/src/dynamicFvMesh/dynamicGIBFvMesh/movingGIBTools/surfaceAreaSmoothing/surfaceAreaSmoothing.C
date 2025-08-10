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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/movingGIBTools/surfaceAreaSmoothing/surfaceAreaSmoothing.H"
#include "faMesh/faPatches/faPatch/faPatchData.H"
#include "faMesh/faPatches/constraint/symmetry/symmetryFaPatch.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "algorithms/PatchEdgeFaceWave/PatchEdgeFaceWave.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "algorithms/PatchEdgeFaceWave/patchEdgeFaceBoolInfo.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

labelList Foam::surfaceAreaSmoothing::labelListPatch(const label& masterID)
{
    labelList labelPatchList (mesh_.boundary()[masterID].size(), 0);
    forAll(mesh_.boundary()[masterID], fI)
    {
        const label& gfI = mesh_.boundary()[masterID].patch().localToGlobal(fI);
        labelPatchList[fI] = gfI;
    }
    return labelPatchList;
}

void Foam::surfaceAreaSmoothing::smooth()
{
    scalar GsmoothedMin;
    scalar GsmoothedMax;
    dimensionedScalar ss("ss", dimless, 1);

    areaScalarField& Gsource = *Gsource_;
    areaScalarField& Gsm = *Gsm_;
    edgeScalarField& DsF = *DsF_;
    volScalarField& Gsmoothed = *Gsmoothed_;
    geodeticWallDist& y = *y_;

    if (smoothTransition_)
    {
        scalar zeroBand = 0.7*r_.value();
        forAll(areaMesh_.faces(), fI)
        {
            if (y.y()[fI]<zeroBand)
            {
                Gsource[fI] = 0;
            }
        }
    }




    scalar initMass = gSum(mag(Gsource.primitiveField()*areaMesh_.S()));

    for (int iS=0; iS<5; iS++)
    {
        faMatrix<scalar> GsmEqn
        (
            - fam::laplacian(DsF, Gsm) + fam::Sp(ss, Gsm) - Gsource
        );
        GsmEqn.solve();
    }

    scalar finMass = gSum(mag(Gsm.primitiveField()*areaMesh_.S()));
    GsmoothedMax = max(Gsm.primitiveField());
    GsmoothedMin = min(Gsm.primitiveField());
    reduce(std::tie(GsmoothedMax, GsmoothedMin), ParallelOp<maxOp<scalar>, minOp<scalar>>{});

    if (true)
    {
        scalar fract = 0;
        if (finMass>0)
        {
            fract = initMass/finMass;
        }
        Gsm *= fract;
        finMass = gSum(mag(Gsm.primitiveField())*areaMesh_.S());
    }

    volSurfaceMapping vsm(areaMesh_);
    Gsmoothed.boundaryFieldRef()[patchID_] = Gsm;

    if (false)
    {
//        volScalarField distFromWall = wallDist(mesh_);
        volScalarField distFromWall  (wallDist::New(this->mesh_).y());
        const labelList& fcs = mesh_.boundary()[patchID_].faceCells();
        volScalarField::Boundary& distFromWallbf =
            distFromWall.boundaryFieldRef();
        forAll(distFromWallbf[patchID_], pfI)
        {
            distFromWallbf[patchID_][pfI] = distFromWall[fcs[pfI]];
        }

    //- only for debugging purposes
        const fvPatch& gibPatch(mesh_.boundary()[patchID_]);
        indirectPolyPatch gibPolyPatch =
                    refCast<const indirectPolyPatch>(gibPatch.patch());
        faceList facesF = gibPatch.patch().localFaces();
        pointField pointsF = gibPatch.patch().localPoints();

        PrimitivePatch<faceList, const pointField&> pp =
            PrimitivePatch<faceList, const pointField&>
            (
                facesF,
                pointsF
            );

        simpleVTKWriter bla = simpleVTKWriter
        (
        pp,
        pp.points()
        );
        bla.addFaceData("yy", y.y());
        bla.addFaceData("Gs", Gsource);
        bla.addFaceData("sens", G_.boundaryField()[patchID_]);
        bla.addFaceData("Gsmooth", Gsmoothed.boundaryField()[patchID_]);
        bla.addFaceData("wallDist", distFromWall.boundaryField()[patchID_]);
        bla.write
        (
            word("Test_yy") +Gsource.mesh().time().timeName()+".vtk"
        );
    }
}


void Foam::surfaceAreaSmoothing::createFAFields()
{

    scalar sigmaFactor = 2*standardDeviation_*standardDeviation_;
    dimensionedScalar Ds(r_*r_/sigmaFactor);

    DsF_ = new edgeScalarField
    (
        IOobject
        (
            "diffusivity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        areaMesh_,
        Ds
    );


    Gsmoothed_ = new volScalarField
    (
        IOobject
        (
            "Gsmoothed",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Gsmooth", dimless, 0.0)
    );

    Gsm_ = new areaScalarField
    (
        IOobject
        (
            "Gsm",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        areaMesh_,
        dimensionedScalar("Gsm", dimless, 0),
        "fixedValue"
    );

    Gsource_ = new areaScalarField
    (
        IOobject
        (
            "Gsource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        areaMesh_,
        dimensionedScalar("Gsource", dimless, 0)
    );

    areaScalarField& Gsource = *Gsource_;

    areaScalarField& Gsm = *Gsm_;

    edgeScalarField& DsF = *DsF_;

    scalarField curv = mesh_.boundary()[patchID_].patch().faceCurvature();
    const labelListList& edges = areaMesh_.patch().faceEdges();

    edgeScalarField::Boundary& DsFbf = DsF.boundaryFieldRef();
    forAll(areaMesh_.faces(), fI)
    {
        if (!constraintFaces_[fI])
        {
            Gsource[fI] = G_.boundaryField()[patchID_][fI];
        }
        else
        {
            const labelList& f = edges[fI];
            Gsource[fI] = G_.boundaryField()[patchID_][fI];

            forAll(f, eI)
            {
                const label edgeI = f[eI];
                if (areaMesh_.isInternalEdge(edgeI))
                {
                    DsF[edgeI] = 0;
                }
                else
                {
                    label patchi = areaMesh_.boundary().whichPatch(edgeI);
                    if (areaMesh_.boundary()[patchi].size())
                    {
                        label e = areaMesh_.boundary()[patchi].whichEdge(edgeI);
                        DsFbf[patchi][e] = 0;
                    }
                }
            }
        }
    }

    Gsm.boundaryFieldRef() = 0;
    dimensionedScalar ss("ss", dimless, 1);

    labelHashSet pIDs;
    forAll(areaMesh_.boundary(), pI)
    {
        if (areaMesh_.boundary()[pI].type() == "patch")
        {
            pIDs.insert(pI);
        }
    }


    y_ = new geodeticWallDist(areaMesh_, pIDs, findIsolatedAreas());
}



void Foam::surfaceAreaSmoothing::createFAMesh()
{
    const indirectPrimitivePatch& patch = areaMesh_.patch();
    const label nTotalEdges = patch.nEdges();
    const label nInternalEdges = patch.nInternalEdges();
    const labelListList& edgeFaces = mesh_.edgeFaces();

    labelList faceLabels =  areaMesh_.faceLabels();

    labelList faceCells(faceLabels.size(), -1);

    forAll(faceCells, faceI)
    {
        label faceID = faceLabels[faceI];

        faceCells[faceI] = mesh_.faceOwner()[faceID];
    }

    labelList meshEdges =
        patch.meshEdges
        (
            mesh_.edges(),
            mesh_.cellEdges(),
            faceCells
        );

    labelList bndEdgeFaPatchIDs(nTotalEdges - nInternalEdges, -1);
    labelList isAProcEdge(bndEdgeFaPatchIDs.size(), 0);
    labelList bndEdgeList(bndEdgeFaPatchIDs.size(), -1);


    enum patchType
    {
        processorPatch = 1,
        fixedValue = -2,
        symmetry = -3,
        empty = -4
    };


    const indirectPolyPatch& gibPolyPatch =
        refCast<const
        indirectPolyPatch>(mesh_.boundary()[patchID_].patch());
    const labelList& addr = gibPolyPatch.fAddr();
    boolList gibFaces (mesh_.nFaces(), false);
    forAll(addr, fI)
    {
        gibFaces[addr[fI]] = true;
    }


    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; edgeI++)
    {
        const label& gEdgeI = meshEdges[edgeI];
        label bedgeI = edgeI-nInternalEdges;

        bndEdgeList[bedgeI] = gEdgeI;

        const labelList& edgeFacesI = edgeFaces[gEdgeI];

        forAll(edgeFacesI, faceI)
        {
            label curFace = edgeFacesI[faceI];
            label curPatchID = mesh_.boundaryMesh().whichPatch(curFace);
            if (curPatchID != -1)
            {
                if (mesh_.boundaryMesh().types()[curPatchID]=="empty")
                {
                    bndEdgeFaPatchIDs[bedgeI] = empty;
                }
            }
        }

        forAll(edgeFacesI, faceI)
        {
            label curFace = edgeFacesI[faceI];
            label curPatchID = mesh_.boundaryMesh().whichPatch(curFace);
            if ((curPatchID != -1) && (bndEdgeFaPatchIDs[bedgeI]==-1))
            {
                if (mesh_.boundaryMesh().types()[curPatchID]=="symmetryPlane"
                    || mesh_.boundaryMesh().types()[curPatchID]=="symmetry")
                {
                    bndEdgeFaPatchIDs[bedgeI] = symmetry;
                }
            }
        }
        forAll(edgeFacesI, faceI)
        {
            label curFace = edgeFacesI[faceI];
            label curPatchID = mesh_.boundaryMesh().whichPatch(curFace);
            if ((curPatchID != -1) && (bndEdgeFaPatchIDs[bedgeI]==-1))
            {
                const polyPatch& pp = mesh_.boundaryMesh()[curPatchID];
                if (pp.coupled())
                {
                    bndEdgeFaPatchIDs[bedgeI] = curPatchID;
                    isAProcEdge[bedgeI] = 1;
                }
            }
        }
        if (bndEdgeFaPatchIDs[bedgeI]==-1)
        {
            bndEdgeFaPatchIDs[bedgeI] = fixedValue;
        }
    }

    //Sync edges that belong to proc patch but in the finite area mesh are boundary
    //edges

    syncTools::syncEdgeList
    (
        mesh_,
        bndEdgeList,
        isAProcEdge,
        minEqOp<label>(),
        label(0)
    );


    forAll(bndEdgeList, edgeI)
    {
        label patch = bndEdgeFaPatchIDs[edgeI];
        if
        (
            patch != fixedValue
            && patch != symmetry
            && patch != empty
            && patch>=0
            && !isAProcEdge[edgeI]
        )
        {
            if (mesh_.boundaryMesh()[patch].coupled())
            {
                bndEdgeFaPatchIDs[edgeI] = -2;
            }
        }
    }

    //Find number of total faPathces
    labelList sL = bndEdgeFaPatchIDs;
       DynamicList<label> faPatchID;
    if (sL.size()>0)
    {
        sort(sL);
        faPatchID.append(sL[0]);
        for (int i=0; i < sL.size()-1; i++)
        {
            if (sL[i]<sL[i+1])
            {
                faPatchID.append(sL[i+1]);
            }
        }
    }
    else
    {
        faPatchID.append(-2);
    }
    List<faPatchData> faPatches(faPatchID.size());

    //build faPatchData

    forAll(faPatches, pI)
    {
        SLList<label> tmpList;

        forAll(bndEdgeFaPatchIDs, eI)
        {
            if (bndEdgeFaPatchIDs[eI]==faPatchID[pI])
            {
                tmpList.append(nInternalEdges+eI);
            }
        }
        faPatches[pI].edgeLabels_ = tmpList;
        if (faPatchID[pI] == fixedValue)
        {
            faPatches[pI].name_ = "constBounds";
            faPatches[pI].type_ = "patch";
        }
        else if (faPatchID[pI] == symmetry)
        {
            faPatches[pI].name_ = "symmetryPlane";
            faPatches[pI].type_ = "symmetry";
        }
        else if (faPatchID[pI] == empty)
        {
            faPatches[pI].name_ = "emptyPatches";
            faPatches[pI].type_ = "empty";
        }
        else
        {
               faPatches[pI].name_ = mesh_.boundaryMesh()[faPatchID[pI]].name();
            faPatches[pI].type_ = processorFaPatch::typeName;
            faPatches[pI].ngbPolyPatchID_ = faPatchID[pI];
        }
    }

    // Reorder processorFaPatch using
    // ordering of ngb processorPolyPatch
    forAll(faPatches, patchI)
    {
        if (faPatches[patchI].type_ == processorFaPatch::typeName)
        {
            labelList ngbFaces(faPatches[patchI].edgeLabels_.size(), -1);

            forAll(ngbFaces, edgeI)
            {
                label curEdge = faPatches[patchI].edgeLabels_[edgeI];

                label curPMeshEdge = meshEdges[curEdge];

                forAll(edgeFaces[curPMeshEdge], faceI)
                {
                    label curFace = edgeFaces[curPMeshEdge][faceI];

                    label curPatchID =
                        mesh_.boundaryMesh().whichPatch(curFace);

                    if (curPatchID == faPatches[patchI].ngbPolyPatchID_)
                    {
                        ngbFaces[edgeI] = curFace;
                    }
                }
            }
            SortableList<label> sortedNgbFaces(ngbFaces);
            labelList reorderedEdgeLabels(ngbFaces.size(), -1);
            for (label i = 0; i < reorderedEdgeLabels.size(); i++)
            {
                reorderedEdgeLabels[i] =
                    faPatches[patchI].edgeLabels_
                    [
                        sortedNgbFaces.indices()[i]
                    ];
            }
            faPatches[patchI].edgeLabels_ = reorderedEdgeLabels;
        }
    }

    List<faPatch*> faPatchLst(faPatchID.size());

    forAll(faPatchLst, pI)
    {
        faPatches[pI].dict_.add("type", faPatches[pI].type_);
        faPatches[pI].dict_.add("edgeLabels", faPatches[pI].edgeLabels_);
        if (faPatches[pI].type_ == processorFaPatch::typeName)
        {
            const processorPolyPatch& procPatch
                = refCast<const processorPolyPatch>(mesh_.boundaryMesh()[faPatchID[pI]]);
            faPatches[pI].dict_.add("myProcNo", procPatch.myProcNo());
            faPatches[pI].dict_.add
            (
                "neighbProcNo",
                procPatch.neighbProcNo()
            );
            faPatches[pI].dict_.add
            (
                "ngbPolyPatchIndex",
                faPatchID[pI]
            );
        }
        else if (faPatches[pI].type_ == symmetryFaPatch::typeName)
        {
            faPatches[pI].dict_.add
            (
                "ngbPolyPatchIndex",
                100001
            );
        }
        else
        {
            faPatches[pI].dict_.add
            (
                "ngbPolyPatchIndex",
                100000
            );
        }

        faPatchLst[pI] =
            faPatch::New
            (
                faPatches[pI].name_,
                faPatches[pI].dict_,
                pI,
                areaMesh_.boundary()
            ).ptr();
    }

    areaMesh_.addFaPatches(faPatchLst);

    Pstream::barrier();
}


labelList Foam::surfaceAreaSmoothing::findIsolatedAreas()
{
    const label nTotalEdges = areaMesh_.nEdges();

    DynamicList<label> initialEdges(nTotalEdges);
    DynamicList<patchEdgeFaceBoolInfo> initialEdgesInfo(nTotalEdges);

    boolList markedEdges(areaMesh_.nEdges(), false);

    forAll(areaMesh_.boundary(), patchI)
    {
        if (areaMesh_.boundary()[patchI].type() == "patch")
        {
            labelList pEdges = areaMesh_.boundary()[patchI];
            forAll(pEdges, eI)
            {
                markedEdges[pEdges[eI]] = true;
            }
        }
    }

    forAll(markedEdges, eI)
    {
        if (markedEdges[eI])
        {
            initialEdges.append(eI);
            initialEdgesInfo.append
            (
                patchEdgeFaceBoolInfo
                (
                    true
                )
            );
        }
    }

    initialEdges.shrink();
    initialEdgesInfo.shrink();

    List<patchEdgeFaceBoolInfo> allEdgeInfo(areaMesh_.nEdges());
    List<patchEdgeFaceBoolInfo> allFaceInfo(areaMesh_.nFaces());

    const indirectPolyPatch& gibPolyPatch =
        refCast<const
        indirectPolyPatch>(mesh_.boundary()[patchID_].patch());

    // Walk
    PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchEdgeFaceBoolInfo
    > calc
    (
        mesh_,
        gibPolyPatch,
        initialEdges,
        initialEdgesInfo,
        allEdgeInfo,
        allFaceInfo,
        returnReduce(areaMesh_.nEdges(), sumOp<label>())
    );

    DynamicList<label> dfl(allFaceInfo.size());
    forAll(allFaceInfo, faceI)
    {
        if (!allFaceInfo[faceI].origin())
        {
            dfl.append(faceI);
        }
    }
    return labelList(dfl, true);
}


void Foam::surfaceAreaSmoothing::setSmoothingOptions
(
    const dictionary& dict
)
{
    r_ = dimensionedScalar
    (
        "radius",
        dimensionSet(0,1,0,0,0,0,0),
        dict.lookup("radius")
    );
    standardDeviation_ = dict.lookupOrDefault("Sigma_deviation", 3);
    cutOffAngle_ = dict.lookupOrDefault<scalar>("cut_off_angle", 180);
    smoothTransition_ = dict.lookupOrDefault<bool>
        ("smoothTransition", true);
}

void Foam::surfaceAreaSmoothing::clearOutData()
{
    deleteDemandDrivenData(Gsm_);
    deleteDemandDrivenData(Gsource_);
    deleteDemandDrivenData(DsF_);
    deleteDemandDrivenData(Gsmoothed_);
    deleteDemandDrivenData(y_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceAreaSmoothing::surfaceAreaSmoothing
(
    const fvMesh& mesh,
    const dictionary& dict,
    const label& masterID,
    boolList constraintFaces
)
:
    mesh_(mesh),
    areaMesh_
    (
        mesh_,
        labelListPatch
        (
            masterID
        )
    ),
    patchID_(masterID),
    constraintFaces_(constraintFaces),
    r_
    (
        "radius",
        dimensionSet(0,1,0,0,0,0,0),
        dict.lookup("radius")
    ),
    standardDeviation_(dict.lookupOrDefault("Sigma_deviation", 3)),
    cutOffAngle_(dict.lookupOrDefault<scalar>("cut_off_angle", 180)),
    smoothTransition_
    (
        dict.lookupOrDefault<bool>("smoothTransition", true)
    ),
    G_(mesh.lookupObject<volScalarField>("G")),
    Gsm_(nullptr),
    Gsource_(nullptr),
    DsF_(nullptr),
    Gsmoothed_(nullptr),
    y_(nullptr)
{
    createFAMesh();
    createFAFields();
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::surfaceAreaSmoothing::update()
{
    smooth();
}

tmp<scalarField> Foam::surfaceAreaSmoothing::smoothSens() const
{
    const areaScalarField& Gsm = *Gsm_;
    return Gsm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
