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
    (c) 2017-2023 Esi Ltd.

Application
    transformPoints

Group
    grpMeshManipulationUtilities

Description
    Transforms the mesh points in the polyMesh directory according to the
    translate, rotate and scale options.

Usage
    Supported transformations:

    -translate vector
        Translates the points by the given vector,

    -rotate (vector vector)
        Rotates the points from the first vector to the second,

     or -yawPitchRoll (yawdegrees pitchdegrees rolldegrees)
     or -rollPitchYaw (rolldegrees pitchdegrees yawdegrees)

    -scale vector
        Scales the points by the given vector.

    The any or all of the three options may be specified and are processed
    in the above order.

    With -rotateFields (in combination with -rotate/yawPitchRoll/rollPitchYaw)
    it will also read & transform vector & tensor fields.

    Note:
      roll: rotation about x
      pitch: rotation about y
      yaw: rotation about z

    Options:
      - \par -pointSet \<name\> \n
        Only transform points in the given point set.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/ReadFields/ReadFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "sets/topoSets/pointSet.H"
#include "fields/Fields/transformField/transformField.H"
#include "fields/GeometricFields/transformGeometricField/transformGeometricField.H"
#include "db/IOstreams/StringStreams/IStringStream.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "regionProperties/regionProperties.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void readAndRotateFields
(
    PtrList<GeoField>& flds,
    const fvMesh& mesh,
    const tensor& T,
    const IOobjectList& objects
)
{
    ReadFields(mesh, objects, flds);
    forAll(flds, i)
    {
        Info<< "Transforming " << flds[i].name() << endl;
        dimensionedTensor dimT("t", flds[i].dimensions(), T);
        transform(flds[i], dimT, flds[i]);
    }
}


void rotateFields(const argList& args, const Time& runTime, const tensor& T)
{
    #include "include/createNamedMesh.H"

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    readAndRotateFields(vsFlds, mesh, T, objects);

    PtrList<volVectorField> vvFlds;
    readAndRotateFields(vvFlds, mesh, T, objects);

    PtrList<volSphericalTensorField> vstFlds;
    readAndRotateFields(vstFlds, mesh, T, objects);

    PtrList<volSymmTensorField> vsymtFlds;
    readAndRotateFields(vsymtFlds, mesh, T, objects);

    PtrList<volTensorField> vtFlds;
    readAndRotateFields(vtFlds, mesh, T, objects);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    readAndRotateFields(ssFlds, mesh, T, objects);

    PtrList<surfaceVectorField> svFlds;
    readAndRotateFields(svFlds, mesh, T, objects);

    PtrList<surfaceSphericalTensorField> sstFlds;
    readAndRotateFields(sstFlds, mesh, T, objects);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    readAndRotateFields(ssymtFlds, mesh, T, objects);

    PtrList<surfaceTensorField> stFlds;
    readAndRotateFields(stFlds, mesh, T, objects);

    mesh.write();
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transform (translate/rotate/scale) mesh points.\n"
        "Note: roll=rotation about x, pitch=rotation about y, "
        "yaw=rotation about z"
    );
    argList::addOption
    (
        "translate",
        "vector",
        "translate by the specified <vector> - eg, '(1 0 0)'"
    );
    argList::addOption
    (
        "rotate",
        "(vectorA vectorB)",
        "transform in terms of a rotation between <vectorA> and <vectorB> "
        "- eg, '( (1 0 0) (0 0 1) )'"
    );
    argList::addOption
    (
        "rollPitchYaw",
        "vector",
        "rotate by '(roll pitch yaw)' in degrees"
    );
    argList::addOption
    (
        "yawPitchRoll",
        "vector",
        "rotate by '(yaw pitch roll)' in degrees"
    );
    argList::addBoolOption
    (
        "rotateFields",
        "read and transform vector and tensor fields too"
    );
    argList::addOption
    (
        "scale",
        "vector",
        "scale by the specified amount - eg, '(0.001 0.001 0.001)' for a "
        "uniform [mm] to [m] scaling"
    );
    argList::addOption
    (
        "pointSet",
        "pointSet",
        "Point set to limit the transformation to"
    );

    #include "include/addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in regionProperties"
    );
    argList::addOption
    (
        "regions",
        "list",
        "List of regions to stretch"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    const bool allRegions = args.optionFound("allRegions");
    wordList regionNames;
    wordList regionDirs;
    if (allRegions)
    {
        Info<< "Reconstructing for all regions in regionProperties" << nl
            << endl;
        regionProperties rp(runTime);
        forAllConstIter(HashTable<wordList>, rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, i)
            {
                if (findIndex(regionNames, regions[i]) == -1)
                {
                    regionNames.append(regions[i]);
                }
            }
        }
        regionDirs = regionNames;
    }
    else
    {
        word regionName;
        if (args.optionFound("regions"))
        {
            regionNames = args.optionRead<wordList>("regions");
            regionDirs = regionNames;
        }
        else if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
            regionDirs = regionNames;
        }
        else
        {
            regionNames = wordList(1, fvMesh::defaultRegion);
            regionDirs = wordList(1, word::null);
        }
    }

    forAll(regionNames, regionI)
    {
        const word& regionName = regionNames[regionI];
        const word& regionDir = regionDirs[regionI];

        fileName meshDir(regionDir/polyMesh::meshSubDir);

        Info<< "\n\nTransforming points for region " << regionName << nl
            << endl;

        pointIOField points
        (
            IOobject
            (
                "points",
                runTime.findInstance(meshDir, "points"),
                regionDir/polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        const bool doRotateFields = args.optionFound("rotateFields");

        word pointSetName = word::null;
        const bool doPointSet =
            args.optionReadIfPresent("pointSet", pointSetName);

        if (doRotateFields && doPointSet)
        {
            FatalErrorInFunction
                << "Rotation of fields across the entire mesh, and limiting "
                << "the transformation of points to a set, cannot be done "
                << "simultaneously." << exit(FatalError);
        }

        // This is not actually stringent enough
        if (args.options().empty())
        {
            FatalErrorInFunction
                << "No options supplied, please use one or more of "
                "-translate, -rotate or -scale options."
                << exit(FatalError);
        }

        pointField origPoints;
        labelList setPointIDs;
        if (doPointSet)
        {
            Info<< "Transforming points in pointSet " << pointSetName
                << nl << endl;

            // Save a copy of the original points
            origPoints = points;

            setPointIDs =
                pointSet
                (
                    IOobject
                    (
                        pointSetName,
                        runTime.findInstance(meshDir/"sets", word::null),
                        polyMesh::meshSubDir/"sets",
                        runTime,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                ).toc();

            points = pointField(UIndirectList<point>(points, setPointIDs));
        }

        vector v;
        if (args.optionReadIfPresent("translate", v))
        {
            Info<< "Translating points by " << v << endl;

            points += v;
        }

        if (args.optionFound("rotate"))
        {
            Pair<vector> n1n2
            (
                args.optionLookup("rotate")()
             );
            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);
            tensor T = rotationTensor(n1n2[0], n1n2[1]);

            Info<< "Rotating points by " << T << endl;

            points = transform(T, points);

            if (doRotateFields)
            {
                rotateFields(args, runTime, T);
            }
        }
        else if (args.optionReadIfPresent("rollPitchYaw", v))
        {
            Info<< "Rotating points by" << nl
                << "    roll  " << v.x() << nl
                << "    pitch " << v.y() << nl
                << "    yaw   " << v.z() << nl;

            // Convert to radians
            v *= pi/180.0;

            quaternion R(quaternion::rotationSequence::XYZ, v);

            Info<< "Rotating points by quaternion " << R << endl;
            points = transform(R, points);

            if (doRotateFields)
            {
                rotateFields(args, runTime, R.R());
            }
        }
        else if (args.optionReadIfPresent("yawPitchRoll", v))
        {
            Info<< "Rotating points by" << nl
                << "    yaw   " << v.x() << nl
                << "    pitch " << v.y() << nl
                << "    roll  " << v.z() << nl;

            // Convert to radians
            v *= pi/180.0;

            scalar yaw = v.x();
            scalar pitch = v.y();
            scalar roll = v.z();

            quaternion R = quaternion(vector(0, 0, 1), yaw);
            R *= quaternion(vector(0, 1, 0), pitch);
            R *= quaternion(vector(1, 0, 0), roll);

            Info<< "Rotating points by quaternion " << R << endl;
            points = transform(R, points);

            if (doRotateFields)
            {
                rotateFields(args, runTime, R.R());
            }
        }

        if (args.optionReadIfPresent("scale", v))
        {
            Info<< "Scaling points by " << v << endl;

            points.replace(vector::X, v.x()*points.component(vector::X));
            points.replace(vector::Y, v.y()*points.component(vector::Y));
            points.replace(vector::Z, v.z()*points.component(vector::Z));
        }

        if (doPointSet)
        {
            UIndirectList<point>(origPoints, setPointIDs) = pointField(points);
            points = origPoints;
        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        Info<< "Writing points into directory " << points.path() << nl << endl;
        points.write();
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
