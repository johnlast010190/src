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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Create mesh, no clear-out\n" << endl;
    polyMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );


    Info<< "Cell centres" << nl << mesh.cellCentres() << endl;
    Info<< "Cell volumes" << nl << mesh.cellVolumes() << endl;
    Info<< "Cell shapes" << nl << mesh.cellShapes() << endl;
    Info<< "Cell face centres" << nl << mesh.faceCentres() << endl;

    // Test construct from cellShapes
    {
        pointField points(mesh.points());
        cellShapeList shapes(mesh.cellShapes());

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        faceListList boundaryFaces(pbm.size());
        forAll(pbm, patchi)
        {
            boundaryFaces[patchi] = pbm[patchi];
        }
        wordList boundaryPatchNames(pbm.names());
        PtrList<dictionary> boundaryDicts(pbm.size());
        forAll(pbm, patchi)
        {
            OStringStream os;
            os << pbm[patchi];
            IStringStream is(os.str());
            boundaryDicts.set(patchi, new dictionary(is));
        }

        word defaultBoundaryPatchName = "defaultFaces";
        word defaultBoundaryPatchType = emptyPolyPatch::typeName;

        polyMesh newMesh
        (
            IOobject
            (
                "newMesh",
                runTime.timeName(),
                runTime,
                Foam::IOobject::NO_READ
            ),
            Xfer<pointField>(points),
            shapes,
            boundaryFaces,
            boundaryPatchNames,
            boundaryDicts,
            defaultBoundaryPatchName,
            defaultBoundaryPatchType
        );

        Info<< "New cell centres" << nl << newMesh.cellCentres() << endl;
        Info<< "New cell volumes" << nl << newMesh.cellVolumes() << endl;
        Info<< "New cell shapes" << nl << newMesh.cellShapes() << endl;
        Info<< "New cell face centres" << nl << newMesh.faceCentres() << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
