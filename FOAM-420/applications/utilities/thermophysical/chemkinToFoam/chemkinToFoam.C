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
    (c) 2011-2017 OpenFOAM Foundation

Application
    chemkinToFoam

Group
    grpSurfaceUtilities

Description
    Converts CHEMKINIII thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "chemistryReaders/chemkinReader/chemkinReader.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/StringStreams/OStringStream.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

#include "cfdTools/general/include/fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
//#include "include/setRootCase.H"


    // Increase the precision of the output for JANAF coefficients
    Ostream::defaultPrecision(10);

    argList::addNote
    (
        "Converts CHEMKINIII thermodynamics and reaction data files into\n"
        "OpenFOAM format."
    );
    argList::noParallel();
    argList::noFunctionObjects();
    argList::validArgs.append("CHEMKINFile");
    argList::validArgs.append("CHEMKINThermodynamicsFile");
    argList::validArgs.append("CHEMKINTransport");
    argList::validArgs.append("FOAMChemistryFile");
    argList::validArgs.append("FOAMThermodynamicsFile");

    argList::addBoolOption
    (
        "newFormat",
        "read Chemkin thermo file in new format"
    );

    argList args(argc, argv);

#include "include/createTime.H"

    word instance = runTime.findInstance(
        polyMesh::meshSubDir, "points", IOobject::READ_IF_PRESENT
    );

    bool meshExists(isFile(runTime.path()/instance/polyMesh::meshSubDir/"faces"));

    fvMesh* meshPtr;

    if (meshExists)
    {
        Info<< "Creating mesh ... " << endl;
        meshPtr =
            new fvMesh
            (
                Foam::IOobject
                    (
                        Foam::fvMesh::defaultRegion,
                        runTime.timeName(),
                        runTime,
                        Foam::IOobject::MUST_READ
                    )
            );
    }
    else
    {
        Info<< "mesh not found" << endl;
        Info<< "Creating createSingleCellMesh on the fly" << endl;
        Info<< "Constructing single cell mesh" << nl << endl;

        labelList owner(6, label(0));
        labelList neighbour(0);

        pointField points(8);
        points[0] = vector(0, 0, 0);
        points[1] = vector(1, 0, 0);
        points[2] = vector(1, 1, 0);
        points[3] = vector(0, 1, 0);
        points[4] = vector(0, 0, 1);
        points[5] = vector(1, 0, 1);
        points[6] = vector(1, 1, 1);
        points[7] = vector(0, 1, 1);
        points[7] = vector(0, 1, 1);

        const cellModel& hexa = *(cellModeller::lookup("hex"));
        faceList faces = hexa.modelFaces();

        meshPtr = new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::READ_IF_PRESENT
            ),
            xferMove<Field<vector>>(points),
            faces.xfer(),
            owner.xfer(),
            neighbour.xfer()
        );

        List<polyPatch*> patches(1);

        patches[0] = new emptyPolyPatch
        (
            "boundary",
            6,
            0,
            0,
            meshPtr->boundaryMesh(),
            emptyPolyPatch::typeName
        );
        meshPtr->addFvPatches(patches);
    }
    const fvMesh& mesh = *meshPtr;

    const bool newFormat = args.optionFound("newFormat");

    speciesTable species;

    chemkinReader cr(species, mesh, args[1], args[3], args[2], newFormat);


    OFstream reactionsFile(args[4]);
    reactionsFile
        << "elements" << cr.elementNames() << token::END_STATEMENT << nl << nl;
    reactionsFile
        << "species" << cr.species() << token::END_STATEMENT << nl << nl;
    cr.reactions().write(reactionsFile);


    // Temporary hack to splice the specie composition data into the thermo file
    // pending complete integration into the thermodynamics structure

    OStringStream os;
    cr.speciesThermo().write(os);
    dictionary thermoDict(IStringStream(os.str())());

    wordList speciesList(thermoDict.toc());

    // Add elements
    forAll(speciesList, si)
    {
        dictionary elementsDict("elements");
        forAll(cr.specieComposition()[speciesList[si]], ei)
        {
            elementsDict.add
            (
                cr.specieComposition()[speciesList[si]][ei].name(),
                cr.specieComposition()[speciesList[si]][ei].nAtoms()
            );
        }

        thermoDict.subDict(speciesList[si]).add("elements", elementsDict);
    }

    thermoDict.write(OFstream(args[5])(), false);


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
