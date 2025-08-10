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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Description

\*---------------------------------------------------------------------------*/

#include "faMesh/faBoundaryMesh/faBoundaryMesh.H"
#include "faMesh/faMesh.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "db/IOstreams/Pstreams/PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faBoundaryMesh, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
faBoundaryMesh::faBoundaryMesh
(
    const IOobject& io,
    const faMesh& mesh
)
:
    faPatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    if (readOpt() == IOobject::MUST_READ)
    {
        faPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchI)
        {
            patches.set
            (
                patchI,
                faPatch::New
                (
                    patchEntries[patchI].keyword(),
                    patchEntries[patchI].dict(),
                    patchI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "faBoundaryMesh::polyBoundaryMesh"
            "(const IOobject&, const faMesh&)"
        );

        close();
    }
}


// Construct given size. Patches will be set later
faBoundaryMesh::faBoundaryMesh
(
    const IOobject& io,
    const faMesh& pm,
    const label size
)
:
    faPatchList(size),
    regIOobject(io),
    mesh_(pm)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate the geometry for the patches (transformation tensors etc.)
void faBoundaryMesh::calcGeometry()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);
    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchI)
        {
            operator[](patchI).initGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchI)
        {
            operator[](patchI).calcGeometry(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();
        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchI = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchI).initGeometry(pBufs);
            }
            else
            {
                operator[](patchI).calcGeometry(pBufs);
            }
        }
    }
}


// Return the mesh reference
const faMesh& faBoundaryMesh::mesh() const
{
    return mesh_;
}


lduInterfacePtrsList faBoundaryMesh::interfaces() const
{
    lduInterfacePtrsList interfaces(size());

    forAll(interfaces, patchi)
    {
        if (isA<lduInterface>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
               &refCast<const lduInterface>(this->operator[](patchi))
            );
        }
    }

    return interfaces;
}


// Return a list of patch types
wordList faBoundaryMesh::types() const
{
    const faPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchI)
    {
        t[patchI] = patches[patchI].type();
    }

    return t;
}


// Return a list of patch names
wordList faBoundaryMesh::names() const
{
    const faPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchI)
    {
        t[patchI] = patches[patchI].name();
    }

    return t;
}


label faBoundaryMesh::findPatchID(const word& patchName) const
{
    const faPatchList& patches = *this;

    forAll(patches, patchI)
    {
        if (patches[patchI].name() == patchName)
        {
            return patchI;
        }
    }

    // Patch not found
    return -1;
}


// Return patch index for a given edge label
label faBoundaryMesh::whichPatch(const label edgeIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (edgeIndex >= mesh().nEdges())
    {
        FatalErrorInFunction
            << "given label greater than the number of edges"
            << abort(FatalError);
    }

    if (edgeIndex < mesh().nInternalEdges())
    {
        return -1;
    }

    forAll(*this, patchI)
    {
        label start = mesh_.patchStarts()[patchI];
        label size = operator[](patchI).faPatch::size();

        if
        (
            edgeIndex >= start
         && edgeIndex < start + size
        )
        {
            return patchI;
        }
    }

    // If not in any of above, it's trouble!
    FatalErrorInFunction
        << "error in patch search algorithm"
        << abort(FatalError);

    return -1;
}

Foam::labelList Foam::faBoundaryMesh::findIndices
(
    const keyType& key,
    const bool usePatchGroups
) const
{
    DynamicList<label> indices;

    if (!key.empty())
    {
        if (key.isPattern())
        {
            indices = findStrings(key, this->names());
        }
        else
        {
            // Literal string. Special version of above to avoid
            // unnecessary memory allocations

            indices.setCapacity(1);
            forAll(*this, i)
            {
                if (key == operator[](i).name())
                {
                    indices.append(i);
                    break;
                }
            }
        }
    }

    return labelList(indices, true);
}

bool faBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalEdges();
    const faBoundaryMesh& bm = *this;

    bool boundaryError = false;

    forAll(bm, patchI)
    {
        if (bm[patchI].start() != nextPatchStart)
        {
            boundaryError = true;

            Info
                << "bool faBoundaryMesh::checkDefinition("
                << "const bool report) const : "
                << "Problem with boundary patch " << patchI
                << ".\nThe patch should start on face no " << nextPatchStart
                << " and the boundary file specifies " << bm[patchI].start()
                << "." << nl << endl;
        }

        nextPatchStart += bm[patchI].faPatch::size();
    }

    if (boundaryError)
    {
        SeriousErrorIn
        (
            "bool faBoundaryMesh::checkDefinition("
            "const bool report) const"
        )   << "This mesh is not valid: boundary definition is in error."
            << endl;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Boundary definition OK." << endl;
        }
    }

    return boundaryError;
}


// Correct faBoundaryMesh after moving points
void faBoundaryMesh::movePoints(const pointField& p)
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            forAll(*this, patchI)
            {
                operator[](patchI).initMovePoints(pBufs, p);
            }

            pBufs.finishedSends();

            forAll(*this, patchI)
            {
                operator[](patchI).movePoints(pBufs, p);
            }
        }
        else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
        {
            const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

            // Dummy.
            pBufs.finishedSends();

            forAll(patchSchedule, patchEvali)
            {
                const label patchI = patchSchedule[patchEvali].patch;

                if (patchSchedule[patchEvali].init)
                {
                    operator[](patchI).initMovePoints(pBufs, p);
                }
                else
                {
                    operator[](patchI).movePoints(pBufs, p);
                }
            }
        }
}


void faBoundaryMesh::updateMesh()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchI)
        {
            operator[](patchI).initUpdateMesh(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchI)
        {
            operator[](patchI).updateMesh(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        forAll(patchSchedule, patchEvali)
        {
            const label patchI = patchSchedule[patchEvali].patch;

            if (patchSchedule[patchEvali].init)
            {
                operator[](patchI).initUpdateMesh(pBufs);
            }
            else
            {
                operator[](patchI).updateMesh(pBufs);
            }
        }
    }
}

// writeData member function required by regIOobject
bool faBoundaryMesh::writeData(Ostream& os) const
{
    const faPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches, patchi)
    {
        os  << indent << patches[patchi].name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << patches[patchi] << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_LIST;

    // Check state of IOstream
    os.check("polyBoundaryMesh::writeData(Ostream& os) const");

    return os.good();
}


Ostream& operator<<(Ostream& os, const faBoundaryMesh& bm)
{
    bm.writeData(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
