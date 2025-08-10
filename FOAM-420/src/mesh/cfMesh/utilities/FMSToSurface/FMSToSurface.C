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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description
    Creates surface patches from surface subsets

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/triSurfaceTools/triSurfaceCopyParts/triSurfaceCopyParts.H"
#include "include/demandDrivenData.H"
#include "db/IOstreams/Fstreams/OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void exportFeatureEdges
(
    const triSurf& origSurf,
    const fileName& edgeFileName
)
{
    OFstream file(edgeFileName);

    const pointField& points = origSurf.points();
    labelList newPointLabel(points.size(), -1);
    label nPoints(0);

    const edgeLongList& featureEdges = origSurf.featureEdges();
    forAll(featureEdges, feI)
    {
        const edge& e = featureEdges[feI];

        if (newPointLabel[e[0]] == -1)
            newPointLabel[e[0]] = nPoints++;
        if (newPointLabel[e[1]] == -1)
            newPointLabel[e[1]] = nPoints++;
    }

    pointField pCopy(nPoints);
    forAll(newPointLabel, pI)
    {
        if (newPointLabel[pI] < 0)
            continue;

        pCopy[newPointLabel[pI]] = points[pI];
    }

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << pCopy.size() << " float\n";
    forAll(pCopy, pI)
    {
        const point& p = pCopy[pI];
        file << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
    }

    file << "\nLINES " << featureEdges.size()
         << ' ' << 3*featureEdges.size() << nl;
    forAll(featureEdges, edgeI)
    {
        const edge& e = featureEdges[edgeI];
        file << "2 " << newPointLabel[e[0]]
             << token::SPACE << newPointLabel[e[1]] << nl;
    }
    file << nl;

    if (!file)
        FatalErrorInFunction
            << "Writting of feature edges failed!" << exit(FatalError);
}

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validOptions.insert("exportSubsets", "");
    argList::validOptions.insert("exportFeatureEdges", "");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    fileName outFileNoExt = outFileName.lessExt();
    fileName outExtension = outFileName.ext();

    Info<< "Out file no ext " << outFileNoExt << endl;
    Info<< "Extension " << outExtension << endl;

    //- read the inout surface
    triSurf origSurf(inFileName);

    //- write the surface in the requated format
    origSurf.writeSurface(outFileName);

    //- export surface subsets as separate surface meshes
    if (args.options().found("exportSubsets"))
    {
        DynList<label> subsetIDs;
        origSurf.facetSubsetIndices(subsetIDs);

        triSurfaceCopyParts copyParts(origSurf);

        forAll(subsetIDs, subsetI)
        {
            //- get the name of the subset
            triSurf copySurf;
            wordList subsetName(1);
            subsetName[0] = origSurf.facetSubsetName(subsetIDs[subsetI]);

            //- create a surface mesh corresponding to the subset
            copyParts.copySurface(subsetName, copySurf);

            //- write the mesh on disk
            fileName fName = outFileNoExt+"_facetSubset_"+subsetName[0];
            fName += '.'+outExtension;

            copySurf.writeSurface(fName);
        }
    }

    if (args.options().found("exportFeatureEdges"))
    {
        fileName fName = outFileNoExt+"_featureEdges";
        fName += ".vtk";
        exportFeatureEdges(origSurf, fName);
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
