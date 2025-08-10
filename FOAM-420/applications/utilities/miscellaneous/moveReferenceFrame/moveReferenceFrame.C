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
    (c) 2023 Esi Ltd.

Application
    moveReferenceFrame

Description
    Application serves for visualisation of the expected reference frame motion.
    It has posibility to attach some geometries to its motion of simply move.

    When motion depends on the computed variables it isn't possible to
    pre-compute motion.

    Example moveReferenceFrame dictionary:

    \verbatim
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      moveReferenceFrameDict;
    }

    // Frame name defined in meshObjects
    allMRF
    {
        writeGeometry true;
        geometries ( "wholeCar.stl" );
        outputFormat stl;
        integerNumbering true;
    }
    rightRearWheel_Frame
    {
        writeGeometry true;
        geometries ( "rightWeel.stl" );
        outputFormat stl;
        integerNumbering true;
    }
    \endverbatim

    Note: If the geometry ins't being written the solver can still run through
    the reference frame motions which can be useful in conjunction with
    framesMonitor that can write data about the frame from each timestep.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fields/volFields/volFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "db/IOobjectList/IOobjectList.H"

#include "surfaceFormats/stl/STLsurfaceFormat.H"
#include "stl/STLReader.H"
#include "stl/STLpoint.H"
#include "fields/Fields/transformField/transformField.H"
#include "referenceFrames/coordinateFrame.H"
#include "triSurface/triSurface.H"

using namespace Foam::fileFormats;

const word returnFileName
(
    const scalar time,
    const word& surfaceType,
    const bool integerNumbering
)
{
    char dataFile[256];
    if (integerNumbering)
    {
        sprintf(dataFile, "_%06d.", label(time));
    }
    else
    {
        sprintf(dataFile, "_%.6f.", time);
    }
    return word(dataFile) + surfaceType;

}


void writeGeometryFile
(
    const Time& runTime,
    const word& frameName,
    const word& outputType,
    const word& geomName,
    const triSurface& geometry,
    const bool integerNumbering
)
{
    const fileName triSurfaceDir = runTime.time().constant()/"triSurface";
    // Create file name for the output time
    word timeFileName = geomName;
    if (!isDir(triSurfaceDir/frameName))
    {
        mkDir(triSurfaceDir/frameName);
    }
    timeFileName.removeExt();
    timeFileName +=
        returnFileName
        (
            integerNumbering ? runTime.timeIndex() : runTime.timeOutputValue(),
            outputType,
            integerNumbering
        );

    // Write the geometry in desired format
    Info<< "Writing " << outputType << " file " << timeFileName << endl;
    geometry.write(triSurfaceDir/frameName/timeFileName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"

    Info<< "\nStarting time loop\n" << endl;

    // Creating dummy one cell mesh
    Info<< "Creating dummy one cell mesh" << endl;
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Xfer<pointField>
        (
            pointField
            (
                vectorList
                ({
                    vector(0, 0, 0),
                    vector(1, 0, 0),
                    vector(1, 1, 0),
                    vector(0, 1, 0),
                    vector(0, 0, 1),
                    vector(1, 0, 1),
                    vector(1, 1, 1),
                    vector(0, 1, 1)
                })
            )
        ),
        Xfer<faceList>
        (
            faceList
            ({
                face({0, 4, 6, 2}),
                face({1, 3, 7, 5}),
                face({2, 6, 7, 3}),
                face({0, 1, 5, 4}),
                face({4, 5, 7, 6}),
                face({0, 2, 3, 1})
            })
        ),
        Xfer<cellList>
        (
            cellList({cell(labelList({0, 1, 2, 3, 4, 5}))})
        ),
        false,
        false
    );

    Info<< "\nStarting time loop\n" << endl;
    IOdictionary moveFrameDict
    (
        IOobject
        (
            "moveReferenceFrameDict",
            mesh.time().caseSystem(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const wordList frameNames(moveFrameDict.toc());
    UPtrList<coordinateFrame> frames(frameNames.size());

    // Number of frames
    const label nFrames = frameNames.size();

    // Variables related to writing geometry
    List<Switch> writeGeometry(nFrames);
    List<Switch> integerNumbering(nFrames);
    wordList outputTypes(nFrames);

    // Allow multiple geometries per frame
    List<fileNameList> triFiles(nFrames);
    List<List<triSurface>> triGeomet(nFrames);

    forAll(frameNames, movei)
    {
        const word& frameName = frameNames[movei];
        frames.set(movei, &coordinateFrame::New(mesh, frameName));
        const dictionary& frameSubDict = moveFrameDict.subDict(frameName);
        writeGeometry[movei] =
            frameSubDict.lookupOrDefault<Switch>("writeGeometry", false);
        integerNumbering[movei] =
            frameSubDict.lookupOrDefault<Switch>("integerNumbering", false);
        if (writeGeometry[movei] && Pstream::master())
        {
            outputTypes[movei] =
                frameSubDict.lookupOrDefault<word>("outputFormat", "stl");
            const wordList geometryNames =
                frameSubDict.lookup<wordList>("geometries");

            triGeomet[movei].setSize(geometryNames.size());
            triFiles[movei].setSize(geometryNames.size());
            forAll(geometryNames, geomi)
            {
                triFiles[movei][geomi] =
                    runTime.time().constant()
                   /"triSurface"
                   /geometryNames[geomi];

                triGeomet[movei][geomi] = triSurface(triFiles[movei][geomi]);
                writeGeometryFile
                (
                    runTime,
                    frameName,
                    outputTypes[movei],
                    geometryNames[geomi],
                    triGeomet[movei][geomi],
                    integerNumbering[movei]
                );
            }
        }
        Info<< endl;
    }

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        // Updating all the coordinate frames that are loaded
        coordinateFrame::updateStates(mesh.thisDb());

        // If write is on rotate surfaces and write them
        forAll(frameNames, movei)
        {
            if
            (
                writeGeometry[movei]
             && runTime.writeTime()
             && Pstream::master()
            )
            {
                const dictionary& frameSubDict =
                    moveFrameDict.subDict(frameNames[movei]);
                const wordList geometryNames =
                    frameSubDict.lookup<wordList>("geometries");
                forAll(geometryNames, geomi)
                {
                    triSurface& geometry = triGeomet[movei][geomi];
                    const coordinateFrame& frame = frames[movei];
                    geometry.movePoints
                    (
                        transformPoints
                        (
                            frame.transformation(),
                            geometry.points()
                        )
                    );
                    writeGeometryFile
                    (
                        runTime,
                        frameNames[movei],
                        outputTypes[movei],
                        geometryNames[geomi],
                        geometry,
                        integerNumbering[movei]
                    );

                    // If time zero motion is defined transform points back
                    if (!frame.isIncrementalMotion())
                    {
                        geometry.movePoints
                        (
                            transformPoints
                            (
                                inv(frame.transformation()),
                                geometry.points()
                            )
                        );
                    }
                }
            }
        }
        Info<< endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
