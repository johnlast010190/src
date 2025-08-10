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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

Description
    Abaqus surface reader.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// extract data for each column of a line -
// slightly adapted to the .inp format
static std::string readINPToken
(
    const string& line,
    const size_t& width,
    size_t& index
)
{
    size_t indexStart, indexEnd;

    indexStart = index;

    indexEnd = line.find(',', indexStart);

    index = indexEnd + 1;

    if (indexEnd == std::string::npos)
    {
        indexEnd = line.size();
        index = indexEnd;
    }

    return line.substr(indexStart, indexEnd - indexStart);
}

//Not sure how this works, but it works fine eventually!
static scalar parseINPCoord(const string& s)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign))());
        scalar exponent = readScalar(IStringStream(s.substr(expSign+1))());

        if (s[expSign] == '-')
        {
            exponent = -exponent;
        }
        return mantissa*pow(10, exponent);
    }
    else
    {
        return readScalar(IStringStream(s)());
    }
}

bool triSurface::readINP(const fileName& fName)
{
    // get the .inp Abaqus file
    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    // Coordinates of Abaqus point
    DynamicList<point> points;
    // Abaqus index of point
    DynamicList<label> indices;
    // Faces in terms of Abaqus point indices
    DynamicList<labelledTri> faces;
    // From face group to patch
    Map<label> groupToPatch;
    label nPatches = 0;
    // Name for face group
    Map<word> groupToName;
    // look for all names
    HashTable<label> groupNames;
    // status read NODE or ELEMENT
    string status ="NULL";
    //element type - S3 TRIA or S4 QUAD
    string eltype ="NOEL";
    // patches label
    label groupId = -1;

    while (is.good())
    {
        size_t linei = 0;
        string line;
        is.getLine(line);

        if (line.empty() || (line[0] == '*'&& line[1] == '*'))
        {
            continue;
        }

        //status NODE or ELEMENT
        if (line.find("*NODE")!=std::string::npos)
        {
            status ="NODE";
            continue;
        }
        else if (line.find("*ELEMENT")!=std::string::npos)
        {
            status ="ELEMENT";

            //skip 1st & 2nd keywords
            readINPToken(line, 8, linei);
            readINPToken(line, 8, linei);
            string patchstring = readINPToken(line, 8, linei);

            if (patchstring.find("ELSET=")!=std::string::npos)
            {
                groupId++;
                std::size_t elsetpos = patchstring.find("ELSET=");
                string patchname = patchstring.substr(elsetpos+6);
                groupNames.insert(string::validate<word>(patchstring),groupId);
                Info<< "Name of the patch "<<groupId<<" is "<<patchname<<endl;
            }

            // find element type now, S3 and S4 supported
            if (line.find("S3,", 14)!=std::string::npos)
            {
                eltype = "S3";
            }
            else if (line.find("S4,", 14)!=std::string::npos)
            {
                eltype = "S4";
            }
            else
            {
                Info<< "Error: unknown element type found "<<nl
                     << "at line: "<<is.lineNumber()<<endl;
            }
            continue;
        }

        if (status == "NODE")
        {
            //start from the beginning of the line (maybe redundant)
            linei = 0;
            //this stuff reads label and coordinates for one line
            label index =
            readLabel(IStringStream(readINPToken(line, 8, linei))());
            scalar x = parseINPCoord(readINPToken(line, 8, linei));
            scalar y = parseINPCoord(readINPToken(line, 8, linei));
            scalar z = parseINPCoord(readINPToken(line, 8, linei));

            //Info<<" Node index "<<index<<" coord  "<<x<<", "<<y<<", "<<z<<endl;

            indices.append(index);
            points.append(point(x, y, z));
        }
        else if (status == "ELEMENT") //S4 case
        {
            linei = 0;
            if (eltype=="S4")
            {
                readLabel(IStringStream(readINPToken(line, 8, linei))());
                label a =
                    readLabel(IStringStream(readINPToken(line, 8, linei))());
                label b =
                    readLabel(IStringStream(readINPToken(line, 8, linei))());
                label c =
                    readLabel(IStringStream(readINPToken(line, 8, linei))());
                label d =
                    readLabel(IStringStream(readINPToken(line, 8, linei))());

                // Convert group into patch
                Map<label>::const_iterator iter = groupToPatch.find(groupId);

                label patchI;
                if (iter == groupToPatch.end())
                {
                    patchI = nPatches++;
                    groupToPatch.insert(groupId, patchI);
                    Info<< "patch "
                         << patchI
                         << " => group "
                         << groupId
                         << endl;
                }
                else
                {
                    patchI = iter();
                }

                faces.append(labelledTri(a, b, c, patchI));
                faces.append(labelledTri(c, d, a, patchI));
            }
            else if (eltype=="S3")
            {
                readLabel(IStringStream(readINPToken(line, 8, linei))());
                label a =
                    readLabel(IStringStream(readINPToken(line, 8, linei))());
                label b =
                    readLabel(IStringStream(readINPToken(line, 8, linei))());
                label c =
                    readLabel(IStringStream(readINPToken(line, 8, linei))());

                // Convert group into patch
                Map<label>::const_iterator iter = groupToPatch.find(groupId);

                label patchI;
                if (iter == groupToPatch.end())
                {
                    patchI = nPatches++;
                    groupToPatch.insert(groupId, patchI);
                    Info<< "patch "
                         << patchI
                         << " => group "
                         << groupId
                         << endl;
                }
                else
                {
                    patchI = iter();
                }
                faces.append(labelledTri(a, b, c, patchI));
            }
        }
    }

    points.shrink();
    indices.shrink();
    faces.shrink();

    Info<< "Read triangles:" << faces.size() << " points:" << points.size()
        << endl;

    {
        // Build inverse mapping (index to point)
        Map<label> indexToPoint(2*indices.size());
        forAll(indices, i)
        {
            indexToPoint.insert(indices[i], i);
        }

        // Relabel faces
        forAll(faces, i)
        {
            labelledTri& f = faces[i];
            f[0] = indexToPoint[f[0]];
            f[1] = indexToPoint[f[1]];
            f[2] = indexToPoint[f[2]];
        }
    }

    // Convert groupToPatch to patchList.
    geometricSurfacePatchList patches(nPatches);
    forAllConstIter(Map<word>, groupToName, iter)
    {
        Map<label>::const_iterator it = groupToPatch.find(iter.key());

        if (it != groupToPatch.end())
        {
            label patchI = groupToPatch[iter.key()];
            word name(iter());

            patches[patchI] = geometricSurfacePatch
            (
                "empty",
                name,
                patchI
            );
        }
        else
        {
            Info<< "Not found in groupToPatch: " << iter.key() << nl
            << " ignoring and assuming zero sized" << endl;
        }
    }

    Info<< "patches:" << patches << endl;

    // Transfer DynamicLists to straight ones.
    pointField allPoints(points.xfer());

    // Create triSurface
    *this = triSurface(faces, patches, allPoints, true);

    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
