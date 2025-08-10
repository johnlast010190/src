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
    Performs point relocations in the mesh (smoothing) in order to
    improve quality measures. It does not make the mesh invalied.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/smoothers/geometry/meshOptimizer/meshOptimizer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.clear();

    argList::validOptions.insert("nLoops", "int");
    argList::validOptions.insert("nIterations", "int");
    argList::validOptions.insert("nSurfaceIterations", "int");
    argList::validOptions.insert("qualityThreshold", "scalar");
    argList::validOptions.insert("constrainedCellsSet", "word");

#include "include/setRootCase.H"
#include "include/createTime.H"

    //- read the settings
    label nIterations(50);
    label nLoops(10);
    label nSurfaceIterations(2);
    scalar qualityThreshold(0.1);

    if (args.options().found("nLoops"))
    {
        nLoops = readLabel(IStringStream(args.options()["nLoops"])());
    }
    else
    {
        Info<< "Default number of loops is 10" << endl;
    }

    if (args.options().found("nIterations"))
    {
        nIterations =
            readLabel(IStringStream(args.options()["nIterations"])());
    }
    else
    {
        Info<< "Default number of iterations is 50" << endl;
    }

    if (args.options().found("nSurfaceIterations"))
    {
        nSurfaceIterations =
            readLabel(IStringStream(args.options()["nSurfaceIterations"])());
    }
    else
    {
        Info<< "Default number of surface iterations is 2" << endl;
    }

    if (args.options().found("qualityThreshold"))
    {
        qualityThreshold =
            readScalar(IStringStream(args.options()["qualityThreshold"])());
    }
    else
    {
        Info<< "Using default quality threshold 0.1" << endl;
    }

    word constrainedCellSet;

    if (args.options().found("constrainedCellsSet"))
    {
        constrainedCellSet = args.options()["constrainedCellsSet"];
    }
    else
    {
        Info<< "No constraints applied on the smoothing procedure" << endl;
    }

    //- load the mesh from disk
    polyMeshGen pmg(runTime);
    pmg.read();

    //- construct the smoother
    meshOptimizer mOpt(pmg);

    if (!constrainedCellSet.empty())
    {
        //- lock cells in constrainedCellSet
        mOpt.lockCellsInSubset(constrainedCellSet);

        //- find boundary faces which shall be locked
        labelLongList lockedBndFaces, selectedCells;

        const label sId = pmg.cellSubsetIndex(constrainedCellSet);
        pmg.cellsInSubset(sId, selectedCells);

        boolList activeCell(pmg.cells().size(), false);
        forAll(selectedCells, i)
            activeCell[selectedCells[i]] = true;
    }

    //- clear geometry information before volume smoothing
    pmg.clearAddressingData();

    //- perform optimisation using the laplace smoother and
    mOpt.optimizeMeshFV
    (
        nLoops,
        nLoops,
        nIterations,
        nSurfaceIterations
    );

    //- perform optimisation of worst quality faces
    mOpt.optimizeMeshFVBestQuality(nLoops, qualityThreshold);

    //- check the mesh again and untangl bad regions if any of them exist
    mOpt.untangleMeshFV(nLoops, nIterations, nSurfaceIterations);

    Info<< "Writing mesh" << endl;
    pmg.write();

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
