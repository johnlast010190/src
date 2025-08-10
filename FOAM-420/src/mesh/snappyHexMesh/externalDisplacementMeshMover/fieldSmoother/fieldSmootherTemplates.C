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
    (c) 2015 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class Type>
void Foam::fieldSmoother::minSmoothField
(
    const label nIter,
    const PackedBoolList& isPatchMasterPoint,
    const PackedBoolList& isPatchMasterEdge,
    const labelList& meshEdges,
    const indirectPrimitivePatch& adaptPatch,
    const scalarField& fieldMinMag,
    Field<Type>& field
) const
{
    const edgeList& edges = adaptPatch.edges();
    const labelList& meshPoints = adaptPatch.meshPoints();
    const List<labelList>& edgeFaces = adaptPatch.edgeFaces();


    scalarField edgeWeights(edges.size());
    scalarField invSumWeight(meshPoints.size());

    boolList isBoundaryEdge(mesh_.nEdges(), false);
    boolList isBoundaryPoint(mesh_.nPoints(), false);

    labelList nExternalEdge(mesh_.nEdges(), 0);

    forAll(edges, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        const labelList& eFaces = edgeFaces[edgeI];
        nExternalEdge[meshEdgeI] = eFaces.size();
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nExternalEdge,
        plusEqOp<label>(),
        label(0)              // null value
    );

    forAll(mesh_.edges(), edgeI)
    {
        edge e = mesh_.edges()[edgeI];
        if (nExternalEdge[edgeI] == 1)
        {
            isBoundaryEdge[edgeI] =  true;
            isBoundaryPoint[e[0]] = true;
            isBoundaryPoint[e[1]] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        true              // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        true              // null value
    );

    meshRefinement::calculateEdgeWeights
    (
        mesh_,
        isPatchMasterEdge,
        isBoundaryEdge,
        isBoundaryPoint,
        meshEdges,
        meshPoints,
        edges,
        edgeWeights,
        invSumWeight
    );

    // Get smoothly varying patch field.
    Info<< typeName << " : Smoothing field ..." << endl;

    for (label iter = 0; iter < nIter; iter++)
    {
        Field<Type> average(adaptPatch.nPoints(), Type(Zero));

        meshRefinement::weightedSum
        (
            mesh_,
            isPatchMasterEdge,
            isBoundaryEdge,
            isBoundaryPoint,
            meshEdges,
            meshPoints,
            edges,
            edgeWeights,
            field,
            average
        );
        average *= invSumWeight;

        // Transfer to field
        forAll(field, pointI)
        {
            if (mag(field[pointI]) < SMALL)
            {
                continue;
            }

            //full smoothing neighbours + point value
            average[pointI] = 0.5*(field[pointI]+average[pointI]);

            if
            (
                mag(average[pointI]) >= 1.5*mag(fieldMinMag[pointI])
            )
            {
                field[pointI] = average[pointI];
            }
        }

        // Do residual calculation every so often.
        if ((iter % 10) == 0)
        {
            scalar resid = meshRefinement::gAverage
            (
                isPatchMasterPoint,
                mag(field-average)()
            );
            Info<< "    Iteration " << iter << "   residual " << resid << endl;
        }
    }
}


// ************************************************************************* //
