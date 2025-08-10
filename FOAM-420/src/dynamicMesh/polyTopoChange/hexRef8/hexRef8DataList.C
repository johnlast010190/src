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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/IOobject/IOobject.H"
#include "containers/Lists/UList/UList.H"

#include "polyTopoChange/hexRef8/hexRef8DataList.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "polyTopoChange/hexRef8/refinementHistory.H"
#include "fvMesh/fvMesh.H"
#include "primitives/ops/flipOp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexRef8DataList::hexRef8DataList
(
    PtrList<fvMesh>& myProcMeshes,
    const labelList& meshToProc,
    const labelListList& meshNoInProc
)
:
    meshToProc_(meshToProc),
    meshNoInProc_(meshNoInProc)
{
    bool hasCellLevel = false;
    bool hasPointLevel = false;
    bool hasLevel0Edge = false;
    bool hasHistory = false;

     cellLevelPtrs_.setSize(myProcMeshes.size());
    pointLevelPtrs_.setSize(myProcMeshes.size());
    level0EdgePtrs_.setSize(myProcMeshes.size());
    refHistoryPtrs_.setSize(myProcMeshes.size());

    word hexRef8Instance;

    forAll(myProcMeshes, i)
    {
        const fvMesh& mesh = myProcMeshes[i];

        if (i == 0)
        {
            hexRef8Instance =
                mesh.time().findInstance
                (
                    "polyMesh",
                    "cellLevel",
                    IOobject::READ_IF_PRESENT
                );
        }

        IOobject io
        (
            "dummy",
            hexRef8Instance,
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        {
            IOobject rio(io);
            rio.rename("cellLevel");
            if (rio.typeHeaderOk<labelIOList>(true))
            {
                Info<< "Reading hexRef8 data : " << rio.name() << endl;
                cellLevelPtrs_.set(i, new labelIOList(rio));
                hasCellLevel = true;
            }
            else if (hasCellLevel)
            {
                rio.readOpt() = IOobject::NO_READ;
                cellLevelPtrs_.set
                (
                    i,
                    new labelIOList(rio, labelList(mesh.nCells(), 0))
                );
            }
        }
        {
            IOobject rio(io);
            rio.rename("pointLevel");
            if (rio.typeHeaderOk<labelIOList>(true))
            {
                Info<< "Reading hexRef8 data : " << rio.name() << endl;
                pointLevelPtrs_.set(i, new labelIOList(rio));
                hasPointLevel = true;
            }
            else if (hasPointLevel)
            {
                rio.readOpt() = IOobject::NO_READ;
                pointLevelPtrs_.set
                (
                    i,
                    new labelIOList(rio, labelList(mesh.nPoints(), 0))
                );
            }
            {
                IOobject rio(io);
                rio.rename("level0Edge");
                if (rio.typeHeaderOk<uniformDimensionedScalarField>(true))
                {
                    Info<< "Reading hexRef8 data : " << rio.name() << endl;
                    level0EdgePtrs_.set(i, new uniformDimensionedScalarField(rio));
                    hasLevel0Edge = true;
                }
                else if (hasLevel0Edge)
                {
                    rio.readOpt() = IOobject::NO_READ;
                    level0EdgePtrs_.set
                    (
                        i,
                        new uniformDimensionedScalarField
                        (
                            rio,
                            dimensionedScalar
                            (
                                "zero",
                                dimLength,
                                level0EdgePtrs_[0].value()
                            )
                        )
                    );
                }
            }
            {
                IOobject rio(io);
                rio.rename("refinementHistory");
                if (rio.typeHeaderOk<refinementHistory>(true))
                {
                    Info<< "Reading hexRef8 data : " << rio.name() << endl;
                    refHistoryPtrs_.set(i, new refinementHistory(rio));
                    hasHistory = true;
                }
                else if (hasHistory)
                {
                    rio.readOpt() = IOobject::NO_READ;
                    refHistoryPtrs_.set
                    (
                        i,
                        new refinementHistory(rio, mesh.nCells(), true)
                    );
                }
            }
        }
    } // forAll myProcMeshes
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hexRef8DataList::~hexRef8DataList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hexRef8DataList::distribute
(
    PtrList<mapDistributePolyMesh>& myProcDist
)
{
    const label nMyMeshes = meshNoInProc_[Pstream::myProcNo()].size();

    labelList           constructSizes(nMyMeshes);
    List<labelListList> subMaps(nMyMeshes);
    List<labelListList> constructMaps(nMyMeshes);
    boolList            subHasFlips(nMyMeshes);
    boolList            constructHasFlips(nMyMeshes);

    if (cellLevelPtrs_.set(0))
    {
        forAll(myProcDist, i)
        {
            constructSizes[i]    = myProcDist[i].cellMap().constructSize();
            subMaps[i]           = myProcDist[i].cellMap().subMap();
            constructMaps[i]     = myProcDist[i].cellMap().constructMap();
            subHasFlips[i]       = myProcDist[i].cellMap().subHasFlip();
            constructHasFlips[i] = myProcDist[i].cellMap().constructHasFlip();
        }

        mapDistributeBase::distribute
        (
            constructSizes,
            subMaps,
            subHasFlips,
            constructMaps,
            constructHasFlips,
            cellLevelPtrs_,
            flipOp(),
            meshToProc_,
            meshNoInProc_
        );
    }
    if (pointLevelPtrs_.set(0))
    {
        forAll(myProcDist, i)
        {
            constructSizes[i]    = myProcDist[i].pointMap().constructSize();
            subMaps[i]           = myProcDist[i].pointMap().subMap();
            constructMaps[i]     = myProcDist[i].pointMap().constructMap();
            subHasFlips[i]       = myProcDist[i].pointMap().subHasFlip();
            constructHasFlips[i] = myProcDist[i].pointMap().constructHasFlip();
        }

        mapDistributeBase::distribute
        (
            constructSizes,
            subMaps,
            subHasFlips,
            constructMaps,
            constructHasFlips,
            pointLevelPtrs_,
            flipOp(),
            meshToProc_,
            meshNoInProc_
        );
    }

    // No need to distribute the level0Edge

    //TODO
    //if (refHistoryPtrs_.set(0) && refHistoryPtrs_(0).active())
    //{
    //    refHistoryPtr_().distribute(map);
    //}
} // distribute


void Foam::hexRef8DataList::reorder
(
    PtrList<mapPolyMesh>& myProcMap
)
{
    if (cellLevelPtrs_.set(0))
    {
        forAll(myProcMap, i)
        {
            inplaceReorder
            (
                myProcMap[i].cellMap(),
                dynamic_cast<labelList&>(cellLevelPtrs_[i])
            );
        }
    }

    if (pointLevelPtrs_.set(0))
    {
        forAll(myProcMap, i)
        {
            inplaceReorder
            (
                myProcMap[i].pointMap(),
                dynamic_cast<labelList&>(pointLevelPtrs_[i])
            );
        }
    }
} // reorder



bool Foam::hexRef8DataList::write() const
{
    const label nMyMeshes = meshNoInProc_[Pstream::myProcNo()].size();

    bool ok = true;
    for (label i = 0; i < nMyMeshes; i++)
    {
        if (cellLevelPtrs_(i))
        {
            ok = ok && cellLevelPtrs_[i].write();
        }
        if (pointLevelPtrs_(i))
        {
            ok = ok && pointLevelPtrs_[i].write();
        }
        if (level0EdgePtrs_(i))
        {
            ok = ok && level0EdgePtrs_[i].write();
        }
        if (refHistoryPtrs_(i))
        {
            ok = ok && refHistoryPtrs_[i].write();
        }
    }
    return ok;
}


// ************************************************************************* //
