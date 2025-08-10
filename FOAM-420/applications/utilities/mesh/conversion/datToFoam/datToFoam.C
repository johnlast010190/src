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
    (c) 2014 blueCAPE Lda

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.

Modifications made:
    - Always open the files in binary mode, because of how things work on
      Windows.

Application
    datToFoam

Group
    grpMeshConversionUtilities

Description
    Reads in a datToFoam mesh file and outputs a points file.  Used in
    conjunction with blockMesh.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "global/unitConversion/unitConversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("dat file");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    #include "include/createTime.H"

    #if defined(WIN64) || defined(WIN32)
    // Use binary mode in case we read binary.
    // Causes windows reading to fail if we don't.
    std::ifstream plot3dFile(args.args()[1].c_str(),
                             ios_base::in|ios_base::binary);
    #else
    std::ifstream plot3dFile(args[1].c_str());
    #endif

    string line;
    std::getline(plot3dFile, line);
    std::getline(plot3dFile, line);

    IStringStream Istring(line);
    word block;
    string zoneName;
    token punctuation;
    label iPoints;
    label jPoints;

    Istring >> block;
    Istring >> block;
    Istring >> zoneName;
    Istring >> punctuation;
    Istring >> block;
    Istring >> iPoints;
    Istring >> block;
    Istring >> jPoints;

    Info<< "Number of vertices in i direction = " << iPoints << endl
        << "Number of vertices in j direction = " << jPoints << endl;

    // We ignore the first layer of points in i and j the biconic meshes
    label nPointsij = (iPoints - 1)*(jPoints - 1);

    pointField points(nPointsij, Zero);

    for (direction comp = 0; comp < 2; comp++)
    {
        label p(0);

        for (label j = 0; j < jPoints; j++)
        {
            for (label i = 0; i < iPoints; i++)
            {
                double coord;
                plot3dFile >> coord;

                // if statement ignores the first layer in i and j
                if (i>0 && j>0)
                {
                    points[p++][comp] = coord;
                }
            }
        }
    }

    // correct error in biconic meshes
    forAll(points, i)
    {
        if (points[i][1] < 1e-07)
        {
            points[i][1] = 0.0;
        }
    }

    pointField pointsWedge(nPointsij*2, Zero);

    fileName pointsFile(runTime.constantPath()/"points.tmp");
    OFstream pFile(pointsFile);

    scalar a(degToRad(0.1));
    tensor rotateZ =
        tensor
        (
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, -::sin(a), ::cos(a)
        );

    forAll(points, i)
    {
        pointsWedge[i] = (rotateZ & points[i]);
        pointsWedge[i+nPointsij] = cmptMultiply
        (
            vector(1.0, 1.0, -1.0),
            pointsWedge[i]
        );
    }

    Info<< "Writing points to: " << nl
         << "    " << pointsFile << endl;
    pFile << pointsWedge;

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
