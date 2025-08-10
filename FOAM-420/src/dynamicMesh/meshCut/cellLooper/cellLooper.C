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

\*---------------------------------------------------------------------------*/

#include "meshCut/cellLooper/cellLooper.H"
#include "meshes/polyMesh/polyMesh.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "meshTools/meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellLooper, 0);
    defineRunTimeSelectionTable(cellLooper, word);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cellLooper> Foam::cellLooper::New
(
    const word& type,
    const polyMesh& mesh
)
{
    const auto ctor = ctorTableLookup("set type", wordConstructorTable_(), type);
    return autoPtr<cellLooper>(ctor(mesh));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::cellLooper::getVertFacesNonEdge
(
    const label celli,
    const label edgeI,
    const label vertI
) const
{
    // Get faces connected to startEdge
    label face0, face1;
    meshTools::getEdgeFaces(mesh(), celli, edgeI, face0, face1);

    const labelList& pFaces = mesh().pointFaces()[vertI];

    labelList vertFaces(pFaces.size());
    label vertFacei = 0;

    forAll(pFaces, pFacei)
    {
        label facei = pFaces[pFacei];

        if
        (
            (facei != face0)
         && (facei != face1)
         && (meshTools::faceOnCell(mesh(), celli, facei))
        )
        {
            vertFaces[vertFacei++] = facei;
        }
    }
    vertFaces.setSize(vertFacei);

    return vertFaces;
}


Foam::label Foam::cellLooper::getFirstVertEdge
(
    const label facei,
    const label vertI
) const
{
    const labelList& fEdges = mesh().faceEdges()[facei];

    forAll(fEdges, fEdgeI)
    {
        label edgeI = fEdges[fEdgeI];

        const edge& e = mesh().edges()[edgeI];

        if ((e.start() == vertI) || (e.end() == vertI))
        {
            return edgeI;
        }
    }

    FatalErrorInFunction
        << "Can not find edge on face " << facei
        << " using vertex " << vertI
        << abort(FatalError);

    return -1;
}


Foam::labelList Foam::cellLooper::getVertEdgesNonFace
(
    const label celli,
    const label facei,
    const label vertI
) const
{
    const labelList& exclEdges = mesh().faceEdges()[facei];

    const labelList& pEdges = mesh().pointEdges()[vertI];

    labelList vertEdges(pEdges.size());
    label vertEdgeI = 0;

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        if
        (
            (findIndex(exclEdges, edgeI) == -1)
         && meshTools::edgeOnCell(mesh(), celli, edgeI)
        )
        {
            vertEdges[vertEdgeI++] = edgeI;
        }
    }

    vertEdges.setSize(vertEdgeI);

    return vertEdges;
}


Foam::label Foam::cellLooper::getMisAlignedEdge
(
    const vector& refDir,
    const label celli
) const
{
    const labelList& cEdges = mesh().cellEdges()[celli];

    label cutEdgeI = -1;
    scalar maxCos = -GREAT;

    forAll(cEdges, cEdgeI)
    {
        label edgeI = cEdges[cEdgeI];

        scalar cosAngle = mag(refDir & meshTools::normEdgeVec(mesh(), edgeI));

        if (cosAngle > maxCos)
        {
            maxCos = cosAngle;

            cutEdgeI = edgeI;
        }
    }

    return cutEdgeI;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellLooper::cellLooper(const polyMesh& mesh)
:
    edgeVertex(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellLooper::~cellLooper()
{}


// ************************************************************************* //
