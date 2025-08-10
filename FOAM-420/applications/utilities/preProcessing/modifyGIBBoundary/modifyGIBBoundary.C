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
    (c) 2018-2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "addGIBBC.H"
#include "readFields.H"
#include "removeZone.H"
#include "readOrCreateDictionary.H"
#include "surfZone/surfZone/surfZoneList.H"
#include "MeshedSurface/MeshedSurfaces.H"
#include "orientFaceZone/orientFaceZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("addGIB", "add GIB boundaries");
    argList::addBoolOption("removeGIB", "remove GIB boundaries");
    argList::addBoolOption("constant", "save mesh at the constant folder");
    argList::addOption
    (
        "name",
        "word",
        "specify name of the GIB patch"
    );

    argList::addOption
    (
        "faceZone",
        "word",
        "specify name of the faceZone attached to the GIB"
    );

    argList::addOption
    (
        "triSurfaceName",
        "word",
        "specify name of the surface that the GIB will be created"
    );
    argList::addOption
    (
        "orientByPatches",
        "wordList",
        "specify name of the patches which are connected to the fluid"
    );
    argList::addOption
    (
        "names",
        "(patch0 .. patchN)",
        "name patterns of GIB patches that we want to remove. Used only in -removeGIB."
    );
    argList::addBoolOption
    (
        "changeFields",
        "changing boundary conditions at GIB boundaries"
    );

    argList::addOption
    (
        "orientByPoint",
        "point",
        "outside point"
    );
    argList::addOption
    (
        "exportGeometry",
        "word",
        "name of geometry created by triangulating the faceZone"
    );

#include "include/addDictOption.H"

#include "include/setRootCase.H"
#include "include/createTime.H"


    if (args.optionFound("addGIB")||args.optionFound("removeGIB"))
    {

        Foam::Info<< "Create mesh for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;

        Foam::fvMesh mesh
        (
            Foam::IOobject
            (
                Foam::fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );
        const polyBoundaryMesh& patches = mesh.boundaryMesh();


        if (args.optionFound("constant"))
        {
            mesh.setInstance
            (
                runTime.time().constant()
            );
        }


        if (args.optionFound("addGIB"))
        {

            word faceZoneName("");
            if (args.optionFound("faceZone"))
            {
                faceZoneName = args["faceZone"];
            }

            IOdictionary cGIBDict = readOrCreateDictionary(runTime, args);

            const dictionary dInit = cGIBDict.subDict("initialization");
            const word initType = dInit.lookup("initializationType");


            PtrList<volScalarField> volScalarFields;
            PtrList<volVectorField> volVectorFields;
            PtrList<volSphericalTensorField> volSphericalTensorFields;
            PtrList<volSymmTensorField> volSymmTensorFields;
            PtrList<volTensorField> volTensorFields;
            if (args.optionFound("changeFields"))
            {
                IOobjectList objects(mesh, runTime.timeName());
                readFields(mesh, objects, volScalarFields);
                readFields(mesh, objects, volVectorFields);
                readFields(mesh, objects, volSphericalTensorFields);
                readFields(mesh, objects, volSymmTensorFields);
                readFields(mesh, objects, volTensorFields);
            }


            //- add boundary
            {
                #include "addGIBBoundary.H"

                if
                (
                    initType=="dynamicGIBFvMesh"
                )
                {

                    Info<< "Creating dynamicGIBMesh for initializing the GIB"
                         <<endl;

                    const word dynamicFvMeshTypeName
                    (
                        dInit.lookup("dynamicFvMesh")
                    );

                    const auto ctor = ctorTableLookup("dynamicGIBFvMesh type", dynamicGIBFvMesh::IOobjectConstructorTable_(), dynamicFvMeshTypeName);
                    IOobject io
                        (
                             dynamicFvMesh::defaultRegion,
                             runTime.timeName(),
                             runTime,
                             IOobject::NO_READ,
                             IOobject::NO_WRITE,
                             false
                         );

                    if (args.optionFound("constant"))
                    {
                        io.instance() = runTime.time().constant();
                    }

                    autoPtr<dynamicGIBFvMesh> dGIBMesht
                    (
                        ctor
                        (
                            io,
                            dInit
                        )
                    );
                    dynamicGIBFvMesh& dGIBMesh = dGIBMesht();

                    dGIBMesh.updateInit(faceZoneName);

                    if (args.optionFound("constant"))
                    {
                        dGIBMesh.setInstance
                        (
                            runTime.time().constant()
                        );
                    }

                    dGIBMesh.write();
                    dGIBMesh.boundaryMesh().write();

                    if (dInit.found("exportGeometry"))
                    {
                        word surfaceName = dInit.lookup("exportGeometry");
                        dGIBMesh.writeGeometry(surfaceName);
                    }

                    rm(dGIBMesh.phi().filePath());

                    if (args.optionFound("changeFields"))
                    {
                        Info<< "Changing boundaryFields:" <<endl;
                        const dictionary bfs = cGIBDict.subDict("boundaryFields");
                        forAllConstIter(IDLList<entry>, bfs, iter)
                        {
                            word name(iter().keyword());
                            dictionary bf = iter().dict();
                            addGIBBC<scalar>(bf, name, mesh, dGIBMesh);
                            addGIBBC<vector>(bf, name, mesh, dGIBMesh);
                            addGIBBC<sphericalTensor>(bf, name, mesh, dGIBMesh);
                            addGIBBC<symmTensor>(bf, name, mesh, dGIBMesh);
                            addGIBBC<tensor>(bf, name, mesh, dGIBMesh);
                        }
                    }
                }
                else //- no initialization
                {
                    if (dInit.found("orientByPoint"))
                    {
                        const point outPoint = dInit.lookup("orientByPoint");
                        orientFaceZone orFaceZone(mesh, faceZoneName);
                        orFaceZone.orientByPoint(outPoint);
                    }
                    if (dInit.found("exportGeometry"))
                    {
                        word surfaceName = dInit.lookup("exportGeometry");

                        const label masterId = mesh.boundary().size()-2;
                        const polyPatch& mpp = mesh.boundary()[masterId].patch();

                        MeshedSurface<face> surf
                        (
                            mpp,
                            true
                        );

                        if (Pstream::master())
                        {

                            const fileName surfaceFileName = surfaceName;

                            fileName globalCasePath
                            (
                                surfaceFileName.isAbsolute()
                              ? surfaceFileName
                              : (
                                    mesh.time().processorCase()
                                  ?mesh.time().rootPath()/mesh.time().globalCaseName()/surfaceFileName
                                  : mesh.time().path()/surfaceFileName
                                )
                            );
                            globalCasePath.clean();

                            if (!exists(globalCasePath.path()))
                            {
                                Info<< "Path to surface does not exist." << endl;
                                Info<< "Creating path: "
                                     << globalCasePath.path()
                                     << endl;

                                mkDir(globalCasePath.path());
                            }

                            Info<< "Writing merged surface to "
                                 << globalCasePath << endl;

                            surf.write(globalCasePath);
                        }

                    }

                    if (args.optionFound("changeFields"))
                    {
                        Info<< "changing boundaryFields:" <<endl;
                        const dictionary bfs = cGIBDict.subDict("boundaryFields");
                        forAllConstIter(IDLList<entry>, bfs, iter)
                        {
                            word name(iter().keyword());
                            dictionary bf = iter().dict();
                            addGIBBC<scalar>(bf, name, mesh, mesh);
                            addGIBBC<vector>(bf, name, mesh, mesh);
                            addGIBBC<sphericalTensor>(bf, name, mesh, mesh);
                            addGIBBC<symmTensor>(bf, name, mesh, mesh);
                            addGIBBC<tensor>(bf, name, mesh, mesh);
                        }
                    }
                }
            }
        }
        else
        {
            //- remove boundary
            #include "removeGIBBoundary.H"
        }
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << "End\n" <<endl;

    return(0);
}


// ************************************************************************* //
