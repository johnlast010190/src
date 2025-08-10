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
    Translates FOAM mesh to AVL's FPMA format

\*---------------------------------------------------------------------------*/

#include "db/objectRegistry/objectRegistry.H"
#include "db/Time/Time.H"
#include "utilities/meshes/polyMeshGen/polyMeshGen.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "primitives/strings/fileName/fileName.H"

#include "utilities/dataConversion/foamToFPMA/fpmaMesh.H"
#include "utilities/dataConversion/foamToFPMA/writeMeshFPMA.H"
#include "utilities/helperFunctions/helperFunctions.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeMeshFPMA(const polyMeshGen& mesh, const word& fName)
{
    const Time& time = mesh.returnTime();

    const word postProcDir = "FPMA";

    fileName postProcPath = time.path()/postProcDir;

    if (!Foam::isDir(postProcPath))
    {
        mkDir(postProcPath);
    }

    // Open the Case file
    const fileName fpmaFileName = fName + ".fpma";

    Info<< "Writting mesh into " << fpmaFileName << endl;

/*    OFstream fpmaGeometryFile
    (
        postProcPath/fpmaFileName,
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
*/

    OFstream fpmaGeometryFile(postProcPath/fpmaFileName);

    // Construct the FIRE mesh
    fpmaMesh Mesh(mesh);
    Mesh.write(fpmaGeometryFile);
}

void createFIRESelections(polyMeshGen& mesh)
{
    if (!Pstream::parRun())
        return;

    const faceListPMG& faces = mesh.faces();
    const PtrList<processorBoundaryPatch>& procBoundaries =
        mesh.procBoundaries();

    //- create face selections from proc patches
    forAll(procBoundaries, patchI)
    {
        word sName = "InterFacesToProc";
        sName += help::scalarToText(procBoundaries[patchI].neiProcNo());
        const label sID = mesh.addFaceSubset(sName);

        label faceI = procBoundaries[patchI].patchStart();
        const label end = faceI + procBoundaries[patchI].patchSize();
        for (;faceI<end;++faceI)
            mesh.addFaceToSubset(sID, faceI);
    }

    //- create cell selections
    DynList<label> subsets;
    mesh.faceSubsetIndices(subsets);
    forAll(subsets, subsetI)
    {
        const word sName = mesh.faceSubsetName(subsets[subsetI]);

        if (sName.substr(0, 10) == "processor_")
        {
            const word newName = "Proc" + sName.substr(10, sName.size()-10);

            labelLongList cellsInSubset;
            mesh.cellsInSubset(subsets[subsetI], cellsInSubset);
            const label subsetID = mesh.addCellSubset(newName);
            forAll(cellsInSubset, i)
                mesh.addCellToSubset(subsetID, cellsInSubset[i]);
        }
    }

    //- creating node selections
    boolList bndVertex(mesh.points().size(), false);
    forAll(mesh.boundaries(), patchI)
    {
        label faceI = mesh.boundaries()[patchI].patchStart();
        const label end = faceI + mesh.boundaries()[patchI].patchSize();
        for (;faceI<end;++faceI)
        {
            const face& f = mesh.faces()[faceI];

            forAll(f, pI)
                bndVertex[f[pI]] = true;
        }
    }

    forAll(procBoundaries, patchI)
    {
        word sName = "InterSurfaceEdgesToProc";
        sName += help::scalarToText(procBoundaries[patchI].neiProcNo());
        const label subsetID = mesh.addPointSubset(sName);

        label faceI = procBoundaries[patchI].patchStart();
        const label end = faceI + procBoundaries[patchI].patchSize();
        for (;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                if (bndVertex[f[pI]])
                    mesh.addPointToSubset(subsetID, f[pI]);
            }
        }
    }
}

}

// ************************************************************************* //
