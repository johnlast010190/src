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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2022-2024 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "fields/Fields/Field/SubField.H"
#include "db/dictionary/entry/entry.H"
#include "db/dictionary/dictionary.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(indirectPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, indirectPolyPatch, indirect);
    addToRunTimeSelectionTable(polyPatch, indirectPolyPatch, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::indirectPolyPatch::createGIBFaceZone
(
    const polyBoundaryMesh& bm,
    const word fzName
)
{
    const label zoneId = bm.mesh().faceZones().findZoneID(fzName);
    if (zoneId!=-1)
    {
        return zoneId;
    }
    else
    {
        Info<< "faceZone was not found:" << endl;
        Info<< "Constructing a new faceZone called " << fzName << endl;
        const polyMesh& pMesh = bm.mesh();
        polyMesh& cpMesh = const_cast<polyMesh&>(pMesh);
        label oldZoneSize = cpMesh.faceZones().size();
        cpMesh.faceZones().setSize(oldZoneSize + 1);
        label index = oldZoneSize;
        cpMesh.faceZones().set
        (
            index,
            new faceZone
            (
                fzName,
                labelUList(),
                boolList(),
                index,
                cpMesh.faceZones()
            )
        );
        return index;
    }
}


Foam::labelList Foam::indirectPolyPatch::calculateAddressing
(
    const word& pType,
    const polyBoundaryMesh& bm,
    const faceZone& fz
)
{
    DynamicList<label> addr(fz.size());
    const boolList& fmo = fz.flipMap();
    const label& iFz = bm.mesh().faceNeighbour().size();
    forAll(fz, fI)
    {
        const label& fzI = fz[fI];
        if (fzI < iFz)
        {
            addr.append(fzI);
        }
        else
        {
            if (!fmo[fI])
            {
                if (pType=="master")
                {
                    addr.append(fzI);
                }
            }
            else
            {
                if (pType=="slave")
                {
                    addr.append(fzI);
                }
            }
        }
    }
    return labelList(addr);
}

Foam::boolList Foam::indirectPolyPatch::calculateFm
(
    const word& pType,
    const polyBoundaryMesh& bm,
    const faceZone& fz
)
{
    DynamicList<bool> fmD(fz.size());
    const boolList& fmo = fz.flipMap();
    const label& iFz = bm.mesh().faceNeighbour().size();
    forAll(fz, fI)
    {
        const label& fzI = fz[fI];
        const bool& fmI = fmo[fI];
        if (fzI < iFz)
        {
            fmD.append(fmI);
        }
        else
        {
            if (!fmo[fI])
            {
                if (pType=="master")
                {
                    fmD.append(fmI);
                }
            }
            else
            {
                if (pType=="slave")
                {
                    fmD.append(fmI);
                }
            }
        }
    }
    fmD.shrink();
    return boolList(fmD, true);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::indirectPolyPatch::movePoints(const pointField& p)
{
    indirectPrimitivePatch::clearGeom();
}


void Foam::indirectPolyPatch::movePoints(PstreamBuffers&, const pointField& p)
{
    indirectPrimitivePatch::clearGeom();
}


void Foam::indirectPolyPatch::updateMesh(PstreamBuffers&)
{
    clearAddressing();
    indirectPrimitivePatch::clearOut();
}


void Foam::indirectPolyPatch::updateGIB()
{
    clearAddressing();
    indirectPrimitivePatch::clearOut();
    fAddr_ =
        calculateAddressing
        (
            indirectPolyPatchType(),
            boundaryMesh(),
            boundaryMesh().mesh().faceZones()[zoneId_]
        );
    fm_ =
        calculateFm
        (
            indirectPolyPatchType(),
            boundaryMesh(),
            boundaryMesh().mesh().faceZones()[zoneId_]
        );
    this->resetAddressing(fAddr_);
}


void Foam::indirectPolyPatch::clearGeom()
{
    indirectPrimitivePatch::clearGeom();
}


void Foam::indirectPolyPatch::rename(const wordList& newNames)
{
    name_ = newNames[index_];
}


void Foam::indirectPolyPatch::reorder(const labelUList& oldToNewIndex)
{
    index_ = oldToNewIndex[index_];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::indirectPolyPatch::indirectPolyPatch
(
    const word& name,
    const label size,
    const label zoneI,
    const word zPType,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, bm.mesh().faces().size(), zoneI, zPType, index, bm, patchType),
    indirectPrimitivePatch
    (
        IndirectList<face>
        (
            bm.mesh().faces(),
            calculateAddressing
            (
                zPType,
                bm,
                bm.mesh().faceZones()[zoneI]
            )
            //bm.mesh().faceZones()[zoneI]
        ),
        bm.mesh().points()
    ),
    nbrPatchName_(word::null),
    fAddr_
    (
        calculateAddressing
        (
            zPType,
            bm,
            bm.mesh().faceZones()[zoneI]
        )
    ),
    fm_
    (
        calculateFm
        (
            zPType,
            bm,
            bm.mesh().faceZones()[zoneI]
        )
    ),
    nbrPatchID_(-1),
    zoneId_(zoneI),
    indirectPolyPatchType_(zPType),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    magFaceAreasPtr_(nullptr)
{}


Foam::indirectPolyPatch::indirectPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, bm.mesh().faces().size(), dict, index, bm, patchType),
    indirectPrimitivePatch
    (
        IndirectList<face>
        (
            bm.mesh().faces(),
            calculateAddressing
            (
//                indirectPolyPatchType(),
                dict.lookup("indirectPolyPatchType"),
                bm,
                bm.mesh().faceZones()
                [
                    createGIBFaceZone(bm, dict.lookup("faceZone"))
                ]
            )
        ),
        bm.mesh().points()
    ),
    nbrPatchName_(dict.lookupOrDefault("neighbourPatch", word::null)),
    fAddr_
    (
        calculateAddressing
        (
            dict.lookup("indirectPolyPatchType"),
            bm,
            bm.mesh().faceZones()
            [
                createGIBFaceZone(bm, dict.lookup("faceZone"))
            ]
        )
    ),
    fm_
    (
        calculateFm
        (
            dict.lookup("indirectPolyPatchType"),
            bm,
            bm.mesh().faceZones()
            [
                createGIBFaceZone(bm, dict.lookup("faceZone"))
            ]
        )
    ),
    nbrPatchID_(-1),
    zoneId_
    (
        createGIBFaceZone(bm, dict.lookup("faceZone"))
    ),
    indirectPolyPatchType_
    (
        dict.lookup("indirectPolyPatchType")
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    magFaceAreasPtr_(nullptr)
{}


Foam::indirectPolyPatch::indirectPolyPatch
(
    const indirectPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    indirectPrimitivePatch
    (
        IndirectList<face>
        (
            bm.mesh().faces(),
            calculateAddressing
            (
                pp.indirectPolyPatchType(),
                bm,
                bm.mesh().faceZones()[pp.zoneId()]
            )
        ),
        bm.mesh().points()
    ),
    nbrPatchName_(pp.nbrPatchName()),
    fAddr_
    (
    ),
    fm_(pp.fm()),
    nbrPatchID_(-1),
    zoneId_(pp.zoneId()),
    indirectPolyPatchType_(pp.indirectPolyPatchType()),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    magFaceAreasPtr_(nullptr)
{
    this->start_ = bm.mesh().faces().size();
}


Foam::indirectPolyPatch::indirectPolyPatch
(
    const indirectPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label& fzI
)
:
    polyPatch(pp, bm),
    indirectPrimitivePatch
    (
        IndirectList<face>
        (
            bm.mesh().faces(),
            calculateAddressing
            (
                pp.indirectPolyPatchType(),
                bm,
                bm.mesh().faceZones()[fzI]
            )
        ),
        bm.mesh().points()
    ),
    nbrPatchName_(pp.nbrPatchName()),
    fAddr_
    (
        calculateAddressing
        (
            pp.indirectPolyPatchType(),
            bm,
            bm.mesh().faceZones()[fzI]
        )
    ),
    fm_
    (
        calculateFm
        (
            pp.indirectPolyPatchType(),
            bm,
            bm.mesh().faceZones()[fzI]
        )
    ),
    nbrPatchID_(-1),
    zoneId_(fzI),
    indirectPolyPatchType_(pp.indirectPolyPatchType()),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    magFaceAreasPtr_(nullptr)
{
    this->start_ = bm.mesh().faces().size();
}



Foam::indirectPolyPatch::indirectPolyPatch(const indirectPolyPatch& p)
:
    polyPatch(p),
    indirectPrimitivePatch(p),
    nbrPatchName_(p.nbrPatchName()),
    fAddr_(p.fAddr()),
    fm_(p.fm()),
    nbrPatchID_(-1),
    zoneId_(p.zoneId_),
    indirectPolyPatchType_(p.indirectPolyPatchType_),
    boundaryMesh_(p.boundaryMesh_),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    magFaceAreasPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::indirectPolyPatch::~indirectPolyPatch()
{
    clearAddressing();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::label Foam::indirectPolyPatch::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(nbrPatchName_);

        if (nbrPatchID_ == -1)
        {
            WarningInFunction
                << "Illegal neighbour patch name " << nbrPatchName_
                << endl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic
        const indirectPolyPatch& nbrPatch = refCast<const indirectPolyPatch>
        (
            this->boundaryMesh()[nbrPatchID_]
        );

        if (nbrPatch.nbrPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << nbrPatchName()
                << endl << " but that in return specifies "
                << nbrPatch.nbrPatchName()
                << endl;
        }
    }
    return nbrPatchID_;
}

const Foam::polyBoundaryMesh& Foam::indirectPolyPatch::boundaryMesh() const
{
    return boundaryMesh_;
}


const Foam::vectorField::subField Foam::indirectPolyPatch::faceCentres() const
{
    if (!faceCentresPtr_)
    {
        faceCentresPtr_ = new vectorField(this->size());
        vectorField& fc = *faceCentresPtr_;

        // get reference to global cell centres
        const vectorField& gfc = boundaryMesh_.mesh().faceCentres();
        const labelList& addr = addressing();

        forAll(addr, facei)
        {
            fc[facei] = gfc[addr[facei]];
        }
    }
    vectorField::subField sfc = vectorField::subField(*faceCentresPtr_, this->size(), 0);

    return sfc;
}


const Foam::vectorField::subField Foam::indirectPolyPatch::faceAreas() const
{
    if (!faceAreasPtr_)
    {
        faceAreasPtr_ = new vectorField(this->size());
        vectorField& fa = *faceAreasPtr_;
        const vectorField& gfa = boundaryMesh_.mesh().faceAreas();
        const boolList& fmap = fm();
        const labelList& addr = addressing();

        forAll(addr, facei)
        {
            //if ( addr[facei] < boundaryMesh().mesh().nInternalFaces())
            if (addr[facei] < boundaryMesh().mesh().faceNeighbour().size())
            {
                if (indirectPolyPatchType_ == "master")
                {
                    if (!fmap[facei])
                    {
                        fa[facei] = gfa[addr[facei]];
                    }
                    else
                    {
                        fa[facei] = -gfa[addr[facei]];
                    }
                }
                else if (indirectPolyPatchType_ == "slave")
                {
                    if (fmap[facei])
                    {
                        fa[facei] = gfa[addr[facei]];
                    }
                    else
                    {
                        fa[facei] = -gfa[addr[facei]];
                    }
                }
                else
                {
                    WarningInFunction
                        << "master or slave are supported"
                        << endl;
                }
            }
            else
            {
                fa[facei] = gfa[addr[facei]];
            }
        }
        vectorField::subField sfa = vectorField::subField(fa, this->size(), 0);
        return sfa;
    }
    vectorField::subField sfaa = vectorField::subField(*faceAreasPtr_, this->size(), 0);

    return sfaa;
}


const Foam::scalarField::subField Foam::indirectPolyPatch::magFaceAreas() const
{
    if (!magFaceAreasPtr_)
    {
        magFaceAreasPtr_ = mag(faceAreas()).ptr();
    }
    return scalarField::subField(*magFaceAreasPtr_, this->size(), 0);
}


// Return the patch face neighbour cell centres
Foam::tmp<Foam::vectorField> Foam::indirectPolyPatch::faceCellCentres() const
{
    if (!faceCentresPtr_)
    {
        faceCentresPtr_ = new vectorField(this->size());
        vectorField& cc = *faceCentresPtr_;

        // get reference to global cell centres
        const vectorField& gcc = boundaryMesh_.mesh().cellCentres();

        const labelUList& faceCells = this->faceCells();

        forAll(faceCells, facei)
        {
            cc[facei] = gcc[faceCells[facei]];
        }
    }
    return *faceCentresPtr_;
}


const Foam::labelUList& Foam::indirectPolyPatch::faceCells() const
{
    if (!faceCellsPtr_)
    {
        const labelList& fo = boundaryMesh().mesh().faceOwner();
        const labelList& fn = boundaryMesh().mesh().faceNeighbour();
        const boolList& fmap = fm();
        const labelList& addr = addressing();

        faceCellsPtr_ = new labelList(this->size());
        labelList& fc = *faceCellsPtr_;
        forAll(fc, facei)
        {
            //if (boundaryMesh().mesh().isInternalFace(addr[facei]))
            if (addr[facei] < boundaryMesh().mesh().faceNeighbour().size())
            {
                if (indirectPolyPatchType_ == "master")
                {
                    if (!fmap[facei])
                    {
                        fc[facei] = fo[addr[facei]];
                    }
                    else
                    {
                        fc[facei] = fn[addr[facei]];
                    }
                }
                else if (indirectPolyPatchType_ == "slave")
                {
                    if (fmap[facei])
                    {
                        fc[facei] = fo[addr[facei]];
                    }
                    else
                    {
                        fc[facei] = fn[addr[facei]];
                    }
                }
                else
                {
                    WarningInFunction
                        << "master or slave are supported"
                        << endl;
                }
            }
            else
            {
                if (indirectPolyPatchType_ == "master")
                {
                    if (!fmap[facei])
                    {
                        fc[facei] = fo[addr[facei]];
                    }
                }
                else if (indirectPolyPatchType_ == "slave")
                {
                    if (fmap[facei])
                    {
                        fc[facei] = fo[addr[facei]];
                    }
                }
            }
        }
    }

    return *faceCellsPtr_;
}


const Foam::labelList& Foam::indirectPolyPatch::meshEdges() const
{
    if (!mePtr_)
    {
        mePtr_ =
            new labelList
            (
                boundaryMesh_.mesh().faceZones()[zoneId_].meshEdges()
            );
    }

    return *mePtr_;
}

void Foam::indirectPolyPatch::GIBAgglomerate
(
    Foam::scalarField& faceWeight
) const
{
    const faceZone& zfl =  boundaryMesh_.mesh().faceZones()[zoneId_];
    forAll(zfl, facei)
    {
        if (zfl[facei] < boundaryMesh().mesh().nInternalFaces())
        {
            faceWeight[zfl[facei]] = 0.0;
        }
    }
}

void Foam::indirectPolyPatch::clearAddressing()
{
    deleteDemandDrivenData(faceCellsPtr_);
    deleteDemandDrivenData(mePtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
    deleteDemandDrivenData(magFaceAreasPtr_);
}


void Foam::indirectPolyPatch::write(Ostream& os) const
{
    os.writeEntry("type", type());
    patchIdentifier::write(os);
    os.writeEntry
    (
        "faceZone",
        boundaryMesh_.mesh().faceZones()[zoneId_].name()
    );
    os.writeEntry("indirectPolyPatchType", indirectPolyPatchType());
    os.writeEntry("nFaces", "0");
    os.writeEntry("startFace", this->start());
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::indirectPolyPatch::operator=(const polyPatch& p)
{
    if (isA<indirectPolyPatch>(p))
    {
        const indirectPolyPatch& indpp =
            refCast<const indirectPolyPatch>(p);
        clearAddressing();

        patchIdentifier::operator=(indpp);
        indirectPrimitivePatch::operator=(indpp);
        fAddr_ = indpp.fAddr();
        fm_ = indpp.fm();
        zoneId_ = indpp.zoneId();
        indirectPolyPatchType_ = indpp.indirectPolyPatchType();
        start_ = indpp.start();
    }
}

const Foam::face& Foam::indirectPolyPatch::operator[](const label i) const
{
    const Foam::face& f = indirectPrimitivePatch::operator[](i);
    return f;
}
Foam::face& Foam::indirectPolyPatch::operator[](const label i)
{
    face& f = indirectPrimitivePatch::operator[](i);
    return f;
}
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const indirectPolyPatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& os, const indirectPolyPatch& p");
    return os;
}


// ************************************************************************* //
