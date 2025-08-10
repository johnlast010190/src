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
    (c) 2010-2012 Esi Ltd.

Description
    Nastran surface reader.

    - Reads NTL file
    \verbatim
        GRID          28        10.20269-.030265-2.358-8
    \endverbatim


\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Skip lines
static void skipLines(IFstream& is, label iter)
{
    for (int i=0;i<iter;i++)
    {
        string line;
        is.getLine(line);
    }
}

// Do weird things to extract number
static scalar parseNTLscalarValue(const string& s)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign-1))());
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


// Read a column of a given width from a fixed-format NTL file
static std::string readNTLtoken
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
        indexEnd = indexStart + width;
        index = indexEnd;
    }

    return line.substr(indexStart, indexEnd - indexStart);
}


bool triSurface::readNTL(const fileName& fName)
{
    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    // coordinates of point
    DynamicList<point> points;
    // Nastran index of point
    DynamicList<label> indices;
    // Faces in terms of Nastran point indices
    DynamicList<labelledTri> faces;
    // From face group to patch
    Map<label> groupToPatch;
    label nPatches = 0;
    // Name for face group
    Map<word> groupToName;
    // look for duplicate names
    HashTable<label> duplicateNames;
    // look for all names
    HashTable<label> groupNames;
    // Name for face group (Case no PSHELL cmd)
    Map<word> groupToNameNoPSHELL;

    // A single warning per unrecognized command
    HashSet<word> unhandledCmd;

    while (is.good())
    {
        size_t linei = 2;
        string line;
        is.getLine(line);

        //Can be more automated, read and process all lines accordingly. Maybe later

        // Reading/Skipping Title Card
        if (line.substr(0, 2) == "25")
        {
            skipLines(is,1);
        }
        // Reading/Skipping Summary Data
        else if (line.substr(0, 2) == "26")
        {
            skipLines(is,1);
        }
        // Reading Node Data
        else if (line.substr(0, 2) == " 1")
        {
            // Get ID
            label index =
                readLabel(IStringStream(readNTLtoken(line,8,linei))());

            is.getLine(line);
            linei = 0;

            // Get coordinates
            scalar x = parseNTLscalarValue(readNTLtoken(line, 16, linei));
            scalar y = parseNTLscalarValue(readNTLtoken(line, 16, linei));
            scalar z = parseNTLscalarValue(readNTLtoken(line, 16, linei));

            indices.append(index);
            points.append(point(x, y, z));

            // Ignore rest
            skipLines(is,1);
        }
        // Reading Element Data
        else if (line.substr(0, 2) == " 2")
        {
            // Get ID
            //label index =
                readLabel(IStringStream(readNTLtoken(line,8,linei))());
            // Get shape
            label shape =
                readLabel(IStringStream(readNTLtoken(line,8,linei))());
            // Get KC
            label KC =
                readLabel(IStringStream(readNTLtoken(line,8,linei))());

            if (shape==3)
            {
                is.getLine(line);
                linei = 0;

                // Get nodes number
                //label NN =
                      readLabel(IStringStream(readNTLtoken(line,8,linei))());
                // Get element config
                //label conf =
                      readLabel(IStringStream(readNTLtoken(line,8,linei))());
                // Get PID
                label PID =
                    readLabel(IStringStream(readNTLtoken(line,8,linei))());

                is.getLine(line);
                linei = 0;
                label a = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                label b = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                label c = readLabel(IStringStream(readNTLtoken(line, 8, linei))());

                // Convert group into patch
                Map<label>::const_iterator iter = groupToPatch.find(PID);

                label patchi;
                if (iter == groupToPatch.end())
                {
                    patchi = nPatches++;
                    groupToPatch.insert(PID, patchi);
                    if (debug)
                    Info<< "patch " << patchi << " => group " << PID << endl;
                }
                else
                {
                    patchi = iter();
                }

                faces.append(labelledTri(a, b, c, patchi));
            }
            else if (shape==4)
            {
                is.getLine(line);
                linei = 0;

                // Get nodes number
                //label NN =
                      readLabel(IStringStream(readNTLtoken(line,8,linei))());
                // Get element config
                //label conf =
                      readLabel(IStringStream(readNTLtoken(line,8,linei))());
                // Get PID
                label PID =
                    readLabel(IStringStream(readNTLtoken(line,8,linei))());

                is.getLine(line);
                linei = 0;
                label a = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                label b = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                label c = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                label d = readLabel(IStringStream(readNTLtoken(line, 8, linei))());

                // Convert group into patch
                Map<label>::const_iterator iter = groupToPatch.find(PID);

                label patchi;
                if (iter == groupToPatch.end())
                {
                    patchi = nPatches++;
                    groupToPatch.insert(PID, patchi);
                    if (debug)
                    Info<< "patch " << patchi << " => group " << PID << endl;
                }
                else
                {
                    patchi = iter();
                }

                faces.append(labelledTri(a, b, c, patchi));
                faces.append(labelledTri(c, d, a, patchi));
            }
            else
            {
                Info<< "Bar shape not implemented" << endl;
                skipLines(is,KC);
            }
        }
        // Reading Material Properties
        else if (line.substr(0, 2) == " 3")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Element Properties
        else if (line.substr(0, 2) == " 4")
        {
            is.getLine(line);
        }
        // Reading Coordinate Frames
        else if (line.substr(0, 2) == " 5")
        {
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Distributed Loads (Data Card number varies)
        else if (line.substr(0, 2) == " 6")
        {
            is.getLine(line);
        }
        // Reading Node Forces (Data Card number varies)
        else if (line.substr(0, 2) == " 7")
        {
            is.getLine(line);
        }
        // Reading Node Displacements
        else if (line.substr(0, 2) == " 8")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Node Temperatures
        else if (line.substr(0, 2) == "10")
        {
            is.getLine(line);
        }
        // Reading Element Temperatures
        else if (line.substr(0, 2) == "11")
        {
            is.getLine(line);
        }
        // Reading MPC Data
        else if (line.substr(0, 2) == "14")
        {
            is.getLine(line);
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Nodal Heat Source
        else if (line.substr(0, 2) == "15")
        {
            is.getLine(line);
        }
        // Reading Distributed Heat Source (Data Card number varies)
        else if (line.substr(0, 2) == "16")
        {
            is.getLine(line);
        }
        // Reading Convection Coefficients
        else if (line.substr(0, 2) == "17")
        {
            is.getLine(line);
            is.getLine(line);
        }
        // Reading Emissivity Values
        else if (line.substr(0, 2) == "18")
        {
            is.getLine(line);
            is.getLine(line);
        }
        else if (line.substr(0, 2) == "98")
        {
            skipLines(is,1);
        }
        // Reached end of file
        else if (line.substr(0, 2) == "99")
        {
            break;
        }
        // Not yet implemented
        else
        {
            continue;
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

    /*
    if (rawPSHELL)
    {
        forAllConstIter(Map<word>, groupToNameNoPSHELL, iter)
        {
            Map<label>::const_iterator it = groupToPatch.find(iter.key());

            if (it != groupToPatch.end())
            {
                label patchi = groupToPatch[iter.key()];
                word name(iter());

                patches[patchi] = geometricSurfacePatch
                (
                    "empty",
                    name,
                    patchi
                );
            }
            else
            {
                Info<< "Not found in groupToPatch: " << iter.key() << nl
                    << " ignoring and assuming zero sized" << endl;
            }
        }
    }
    else
    {
    */
        forAllConstIter(Map<word>, groupToName, iter)
        {
            Map<label>::const_iterator it = groupToPatch.find(iter.key());

            if (it != groupToPatch.end())
            {
                label patchi = groupToPatch[iter.key()];
                word name(iter());

                // if name already exists set to pid label
                if (duplicateNames.found(iter()))
                {
                    name = word("pid") + Foam::name(iter.key());
                }

                patches[patchi] = geometricSurfacePatch
                (
                    "empty",
                    name,
                    patchi
                );
            }
            else
            {
                Info<< "Not found in groupToPatch: " << iter.key() << nl
                    << " ignoring and assuming zero sized" << endl;
            }
        }
    //}


    if (debug)
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
