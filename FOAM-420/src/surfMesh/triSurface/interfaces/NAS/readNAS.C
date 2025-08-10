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

    - Uses the Ansa "$ANSA_NAME" or the Hypermesh "$HMNAME COMP" extensions
      to obtain patch names.
    - Handles Nastran short, long, and comma-separated free formats.
    - Properly handles the Nastran compact floating point notation: \n
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

// Do weird things to extract number
static scalar parseNASCoord(const string& s)
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


// Read a column of a given width from either a fixed-format NAS file, or a
// comma-separated free-format NAS file
static std::string readNASToken
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


bool triSurface::readNAS(const fileName& fName)
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

    bool rawPSHELL = true;
    // Ansa tags. Denoted by $ANSA_NAME. These will appear just before the
    // first use of a type. We read them and store the pshell types which
    // are used to name the patches.
    label ansaId = -1;
    word ansaType;
    string ansaName;

    // A single warning per unrecognized command
    HashSet<word> unhandledCmd;

    while (is.good())
    {
        size_t linei = 0;
        string line;
        is.getLine(line);

        // Ansa extension
        if (line.substr(0, 10) == "$ANSA_NAME")
        {
            string::size_type sem0 = line.find (';', 0);
            string::size_type sem1 = line.find (';', sem0+1);
            string::size_type sem2 = line.find (';', sem1+1);

            if
            (
                sem0 != string::npos
             && sem1 != string::npos
             && sem2 != string::npos
            )
            {
                ansaId = readLabel
                (
                    IStringStream(line.substr(sem0+1, sem1-sem0-1))()
                );
                ansaType = line.substr(sem1+1, sem2-sem1-1);


                if (line.substr(0, 18) == "$ANSA_NAME_COMMENT")
                {
                    string::size_type sem3 = line.find (';', sem2+1);
                    if (sem3 != string::npos)
                    {
                        ansaName = line.substr(sem2+1, sem3-sem2-1);
                    }
                    else
                    {
                        ansaName = line.substr(sem2+1, line.size() - sem2 - 1);
                        while (true)
                        {
                            string buf;
                            is.getLine(buf);
                            string::size_type sem4 = buf.find (';', 0);
                            if (sem4 != string::npos)
                            {
                                ansaName += buf.substr(1, sem4 - 1);
                                break;
                            }
                            else
                            {
                                ansaName += buf.substr(0, buf.size());
                            }
                        }
                    }

                    if (ansaType == "PSHELL")
                    {
                        rawPSHELL = false;
                        if (debug)
                        {
                            Pout<< "Found name " << ansaName << " for group "
                                << ansaId << endl;
                        }

                        if
                        (
                            !groupNames.insert
                             (
                                 string::validate<word>(ansaName),
                                 ansaId
                              )
                         )
                        {
                            duplicateNames.insert
                            (
                                string::validate<word>(ansaName),
                                ansaId
                            );
                        }

                        groupToName.insert
                        (
                            ansaId,
                            string::validate<word>(ansaName)
                        );
                    }
                }
                else
                {
                    string nameString;
                    is.getLine(ansaName);
                    if (ansaName[ansaName.size()-1] == '\r')
                    {
                        ansaName = ansaName.substr(1, ansaName.size()-2);
                    }
                    else
                    {
                        ansaName = ansaName.substr(1, ansaName.size()-1);
                    }

                    groupToNameNoPSHELL.insert
                    (
                        ansaId,
                        string::validate<word>(ansaName)
                     );
                }

                // Info<< "ANSA tag for NastranID:" << ansaId
                //     << " of type " << ansaType
                //     << " name " << ansaName << endl;
            }
        }


        // Hypermesh extension
        // $HMNAME COMP                   1"partName"
        if
        (
            line.substr(0, 12) == "$HMNAME COMP"
         && line.find ('"') != string::npos
        )
        {
            label groupId = readLabel
            (
                IStringStream(line.substr(16, 16))()
            );

            IStringStream lineStream(line.substr(32));

            string rawName;
            lineStream >> rawName;

            if (!groupNames.insert(string::validate<word>(rawName),groupId))
            {
                duplicateNames.insert(string::validate<word>(rawName),groupId);
            }

            groupToName.insert(groupId, string::validate<word>(rawName));
            Info<< "group " << groupId << " => " << rawName << endl;
        }


        if (line.empty() || line[0] == '$')
        {
            // Skip empty or comment
            continue;
        }

        // Check if character 72 is continuation
        if (line.size() > 72 && line[72] == '+')
        {
            line = line.substr(0, 72);

            while (true)
            {
                string buf;
                is.getLine(buf);

                if (buf.size() > 72 && buf[72]=='+')
                {
                    line += buf.substr(8, 64);
                }
                else
                {
                    line += buf.substr(8, buf.size()-8);
                    break;
                }
            }
        }

        // Read first word
        word cmd(IStringStream(readNASToken(line, 8, linei))());

        if (cmd == "CTRIA3")
        {
            readNASToken(line, 8, linei);
            label groupId =
                readLabel(IStringStream(readNASToken(line, 8, linei))());
            label a = readLabel(IStringStream(readNASToken(line, 8, linei))());
            label b = readLabel(IStringStream(readNASToken(line, 8, linei))());
            label c = readLabel(IStringStream(readNASToken(line, 8, linei))());

            // Convert group into patch
            Map<label>::const_iterator iter = groupToPatch.find(groupId);

            label patchi;
            if (iter == groupToPatch.end())
            {
                patchi = nPatches++;
                groupToPatch.insert(groupId, patchi);
                Info<< "patch " << patchi << " => group " << groupId << endl;
            }
            else
            {
                patchi = iter();
            }

            faces.append(labelledTri(a, b, c, patchi));
        }
        else if (cmd == "CTRIA3*")
        {
            readNASToken(line, 16, linei);
            label groupId =
                readLabel(IStringStream(readNASToken(line, 16, linei))());
            label a = readLabel(IStringStream(readNASToken(line, 16, linei))());
            label b = readLabel(IStringStream(readNASToken(line, 16, linei))());

            linei = 0;
            is.getLine(line);
            readNASToken(line, 8, linei);
            label c = readLabel(IStringStream(readNASToken(line, 16, linei))());

            // Convert group into patch
            Map<label>::const_iterator iter = groupToPatch.find(groupId);

            label patchi;
            if (iter == groupToPatch.end())
            {
                patchi = nPatches++;
                groupToPatch.insert(groupId, patchi);
                Info<< "patch " << patchi << " => group " << groupId << endl;
            }
            else
            {
                patchi = iter();
            }

            faces.append(labelledTri(a, b, c, patchi));
        }
        else if (cmd == "CQUAD4")
        {
            readNASToken(line, 8, linei);
            label groupId =
                readLabel(IStringStream(readNASToken(line, 8, linei))());
            label a = readLabel(IStringStream(readNASToken(line, 8, linei))());
            label b = readLabel(IStringStream(readNASToken(line, 8, linei))());
            label c = readLabel(IStringStream(readNASToken(line, 8, linei))());
            label d = readLabel(IStringStream(readNASToken(line, 8, linei))());

            // Convert group into patch
            Map<label>::const_iterator iter = groupToPatch.find(groupId);

            label patchi;
            if (iter == groupToPatch.end())
            {
                patchi = nPatches++;
                groupToPatch.insert(groupId, patchi);
                Info<< "patch " << patchi << " => group " << groupId << endl;
            }
            else
            {
                patchi = iter();
            }

            faces.append(labelledTri(a, b, c, patchi));
            faces.append(labelledTri(c, d, a, patchi));
        }
        else if (cmd == "CQUAD4*")
        {
            readNASToken(line, 16, linei);
            label groupId =
                readLabel(IStringStream(readNASToken(line, 16, linei))());
            label a = readLabel(IStringStream(readNASToken(line, 16, linei))());
            label b = readLabel(IStringStream(readNASToken(line, 16, linei))());

            linei = 0;
            is.getLine(line);
            readNASToken(line, 8, linei);
            label c = readLabel(IStringStream(readNASToken(line, 16, linei))());
            label d = readLabel(IStringStream(readNASToken(line, 16, linei))());

            // Convert group into patch
            Map<label>::const_iterator iter = groupToPatch.find(groupId);

            label patchi;
            if (iter == groupToPatch.end())
            {
                patchi = nPatches++;
                groupToPatch.insert(groupId, patchi);
                Info<< "patch " << patchi << " => group " << groupId << endl;
            }
            else
            {
                patchi = iter();
            }

            faces.append(labelledTri(a, b, c, patchi));
            faces.append(labelledTri(c, d, a, patchi));
        }
        else if (cmd == "PSHELL")
        {
            rawPSHELL = false;
            // Read shell type since group gives patchnames
            label groupId =
                readLabel(IStringStream(readNASToken(line, 8, linei))());
            if (groupId == ansaId && ansaType == "PSHELL")
            {
                groupToName.insert(groupId, string::validate<word>(ansaName));
                Info<< "group " << groupId << " => " << ansaName << endl;
            }
        }
        else if (cmd == "PSHELL*")
        {
            rawPSHELL = false;
            // Read shell type since group gives patchnames
            label groupId =
                readLabel(IStringStream(readNASToken(line, 16, linei))());
            if (groupId == ansaId && ansaType == "PSHELL")
            {
                groupToName.insert(groupId, string::validate<word>(ansaName));
                Info<< "group " << groupId << " => " << ansaName << endl;
            }

            linei = 0;
            is.getLine(line);
        }
        else if (cmd == "GRID")
        {
            label index =
                readLabel(IStringStream(readNASToken(line, 8, linei))());
            readNASToken(line, 8, linei);
            scalar x = parseNASCoord(readNASToken(line, 8, linei));
            scalar y = parseNASCoord(readNASToken(line, 8, linei));
            scalar z = parseNASCoord(readNASToken(line, 8, linei));

            indices.append(index);
            points.append(point(x, y, z));
        }
        else if (cmd == "GRID*")
        {
            // Long format is on two lines with '*' continuation symbol
            // on start of second line.
            // Typical line (spaces compacted)
            // GRID*      126   0 -5.55999875E+02 -5.68730474E+02
            // *         2.14897901E+02
            label index =
                readLabel(IStringStream(readNASToken(line, 16, linei))());
            readNASToken(line, 16, linei);
            scalar x = parseNASCoord(readNASToken(line, 16, linei));
            scalar y = parseNASCoord(readNASToken(line, 16, linei));

            linei = 0;
            is.getLine(line);
            if (line[0] != '*')
            {
                FatalErrorInFunction
                    << "Expected continuation symbol '*' when reading GRID*"
                    << " (double precision coordinate) output" << nl
                    << "Read:" << line << nl
                    << "File:" << is.name()
                    << " line:" << is.lineNumber()
                    << exit(FatalError);
            }
            readNASToken(line, 8, linei);
            scalar z = parseNASCoord(readNASToken(line, 16, linei));

            indices.append(index);
            points.append(point(x, y, z));
        }
        else if (unhandledCmd.insert(cmd))
        {
            Info<< "Unhandled Nastran command " << line << nl
                << "File:" << is.name() << " line:" << is.lineNumber() << endl;
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

    forAllConstIter(Map<label>, groupToPatch, iter)
    {
        label patchi = iter.object();
        patches[patchi] = geometricSurfacePatch
            (
                "empty",
                word("pid") + Foam::name(iter.key()),
                patchi
            );
    }

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
