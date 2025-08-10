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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/meshShapes/cellMatcher/cellMatcher.H"

#include "meshes/primitiveMesh/primitiveMesh.H"
#include "containers/HashTables/Map/Map.H"
#include "meshes/meshShapes/face/faceList.H"
#include "primitives/ints/lists/labelList.H"
#include "containers/Lists/ListOps/ListOps.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellMatcher::cellMatcher
(
    const label vertPerCell,
    const label facePerCell,
    const label maxVertPerFace,
    const word& cellModelName
)
:
    localPoint_(100),
    localFaces_(facePerCell),
    faceSize_(facePerCell, -1),
    pointMap_(vertPerCell),
    faceMap_(facePerCell),
    edgeFaces_(2*vertPerCell*vertPerCell),
    pointFaceIndex_(vertPerCell),
    vertLabels_(vertPerCell),
    faceLabels_(facePerCell),
    cellModelName_(cellModelName),
    cellModelPtr_(nullptr)
{
    forAll(localFaces_, facei)
    {
        face& f = localFaces_[facei];

        f.setSize(maxVertPerFace);
    }

    forAll(pointFaceIndex_, vertI)
    {
        pointFaceIndex_[vertI].setSize(facePerCell);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// look for degenerate faces
bool Foam::cellMatcher::findDegenerateFaces
(
    const faceList& faces,
    const labelList& myFaces,
    const labelList& owner,
    const label& numVert,
    const label& cellI
) const
{

    forAll(myFaces, i)
    {
        const face& face = localFaces_[i];
        label prev = 0;
        label next = 0;
        labelHashSet edgeFaces(myFaces.size());

        for
        (
            label fp = 0;
            fp < faceSize_[i];
            fp++
        )
        {
            next =
                nextVert
                (
                    prev,
                    faceSize_[i],
                    !(owner[faceMap_[i]] == cellI)
                );

            label otherFaceI =
                otherFace
                (
                    numVert,
                    face[prev],
                    face[next],
                    i
                );

            if (!edgeFaces.insert(otherFaceI))
            {
                return true;
            }
            prev = next;
        }
    }

    return false;
}


// Create localFaces_ , pointMap_ , faceMap_
Foam::label Foam::cellMatcher::calcLocalFaces
(
    const faceList& faces,
    const labelList& myFaces
)
{
    // Clear map from global to cell numbering
    localPoint_.clear();

    // Renumber face vertices and insert directly into localFaces_
    label newVertI = 0;
    forAll(myFaces, myFacei)
    {
        label facei = myFaces[myFacei];

        const face& f = faces[facei];
        face& localFace = localFaces_[myFacei];

        // Size of localFace
        faceSize_[myFacei] = f.size();

        forAll(f, localVertI)
        {
            label vertI = f[localVertI];

            Map<label>::iterator iter = localPoint_.find(vertI);
            if (iter == localPoint_.end())
            {
                // Not found. Assign local vertex number.

                if (newVertI >= pointMap_.size())
                {
                    // Illegal face: more unique vertices than vertPerCell
                    return -1;
                }

                localFace[localVertI] = newVertI;
                localPoint_.insert(vertI, newVertI);
                newVertI++;
            }
            else
            {
                // Reuse local vertex number.
                localFace[localVertI] = *iter;
            }
        }

        // Create face from localvertex labels
        faceMap_[myFacei] = facei;
    }

    // Create local to global vertex mapping
    forAllConstIter(Map<label>, localPoint_, iter)
    {
        const label fp = iter();
        pointMap_[fp] = iter.key();
    }

    ////debug
    //write(Info);

    return newVertI;
}


void Foam::cellMatcher::calcEdgeAddressing(const label numVert)
{
    edgeFaces_ = -1;

    forAll(localFaces_, localFacei)
    {
        const face& f = localFaces_[localFacei];

        label prevVertI = faceSize_[localFacei] - 1;
        //forAll(f, fp)
        for
        (
            label fp = 0;
            fp < faceSize_[localFacei];
            fp++
        )
        {
            label start = f[prevVertI];
            label end = f[fp];

            label key1 = edgeKey(numVert, start, end);
            label key2 = edgeKey(numVert, end, start);

            if (edgeFaces_[key1] == -1)
            {
                // Entry key1 unoccupied. Store both permutations.
                edgeFaces_[key1] = localFacei;
                edgeFaces_[key2] = localFacei;
            }
            else if (edgeFaces_[key1+1] == -1)
            {
                // Entry key1+1 unoccupied
                edgeFaces_[key1+1] = localFacei;
                edgeFaces_[key2+1] = localFacei;
            }
            else
            {
                FatalErrorInFunction
                    << "edgeFaces_ full at entry:" << key1
                    << " for edge " << start << " " << end
                    << abort(FatalError);
            }

            prevVertI = fp;
        }
    }
}


void Foam::cellMatcher::calcPointFaceIndex()
{
    // Fill pointFaceIndex_ with -1
    forAll(pointFaceIndex_, i)
    {
        labelList& faceIndices = pointFaceIndex_[i];

        faceIndices = -1;
    }

    forAll(localFaces_, localFacei)
    {
        const face& f = localFaces_[localFacei];

        for
        (
            label fp = 0;
            fp < faceSize_[localFacei];
            fp++
        )
        {
            label vert = f[fp];
            pointFaceIndex_[vert][localFacei] = fp;
        }
    }
}


Foam::label Foam::cellMatcher::otherFace
(
    const label numVert,
    const label v0,
    const label v1,
    const label localFacei
) const
{
    label key = edgeKey(numVert, v0, v1);

    if (edgeFaces_[key] == localFacei)
    {
        return edgeFaces_[key+1];
    }
    else if (edgeFaces_[key+1] == localFacei)
    {
        return edgeFaces_[key];
    }
    else
    {
        FatalErrorInFunction
            << "edgeFaces_ does not contain:" << localFacei
            << " for edge " << v0 << " " << v1 << " at key " << key
            << " edgeFaces_[key, key+1]:" <<  edgeFaces_[key]
            << " , " << edgeFaces_[key+1]
            << abort(FatalError);

        return -1;
    }
}


void Foam::cellMatcher::write(Foam::Ostream& os) const
{
    os  << "Faces:" << endl;

    forAll(localFaces_, facei)
    {
        os  << "    ";

        for (label fp = 0; fp < faceSize_[facei]; fp++)
        {
            os  << ' ' << localFaces_[facei][fp];
        }
        os  << endl;
    }

    os  <<  "Face map  : " << faceMap_ << endl;
    os  <<  "Point map : " << pointMap_ << endl;
}


// ************************************************************************* //
