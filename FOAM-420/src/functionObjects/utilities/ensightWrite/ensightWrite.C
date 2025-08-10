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
    (c) 2016-2020 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "ensightWrite/ensightWrite.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "primitives/strings/wordRes/wordRes.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ensightWrite, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ensightWrite,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int Foam::functionObjects::ensightWrite::process(const word& fieldName)
{
    int state = 0;

    writeVolField<scalar>(meshSubset_,fieldName, state);
    writeVolField<vector>(meshSubset_,fieldName, state);
    writeVolField<sphericalTensor>(meshSubset_,fieldName, state);
    writeVolField<symmTensor>(meshSubset_,fieldName, state);
    writeVolField<tensor>(meshSubset_,fieldName, state);

    return state;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ensightWrite::ensightWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    meshSubset_(mesh_),
    writeOpts_
    (
        dict.found("format")
      ? IOstream::formatEnum(dict.lookup("format"))
      : runTime.writeFormat()
    ),
    caseOpts_(writeOpts_.format()),
    selectFields_(),
    selection_(),
    dirName_("ensightWrite"),
    consecutive_(false),
    meshState_(polyMesh::TOPO_CHANGE),
    ensCase_(nullptr),
    ensMesh_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ensightWrite::~ensightWrite()
{
    if (ensCase_.valid())
    {
        // finalize case
        ensCase().write();
        ensCase_.clear();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::ensightWrite::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    readSelection(dict);

    //
    // writer options
    //

    writeOpts_.useBoundaryMesh(dict.lookupOrDefault("boundary", true));
    writeOpts_.useInternalMesh(dict.lookupOrDefault("internal", true));

    // Warn if noPatches keyword (1806) exists and contradicts our settings
    // Cannot readily use Compat since the boolean has the opposite value.
    if
    (
        dict.lookupOrDefault("noPatches", false)
     && writeOpts_.useBoundaryMesh()
    )
    {
        WarningInFunction
            << "Use 'boundary' instead of 'noPatches' to enable/disable "
            << "conversion of the boundaries" << endl;
    }


    if (dict.found("patches"))
    {
        wordReList lst(dict.lookup("patches"));
        wordRes::inplaceUniq(lst);

        writeOpts_.patchSelection(lst);
    }

    if (dict.found("faceZones"))
    {
        wordReList lst(dict.lookup("faceZones"));
        wordRes::inplaceUniq(lst);

        writeOpts_.faceZoneSelection(lst);
    }

    //
    // case options
    //
    caseOpts_.width(dict.lookupOrDefault<label>("width", 8));

    // remove existing output directory
    caseOpts_.overwrite(dict.lookupOrDefault<Switch>("overwrite", false));

    //
    // other options
    //
    dict.readIfPresent("directory", dirName_);
    consecutive_ = dict.lookupOrDefault<Switch>("consecutive", false);

    //
    // output fields
    //
    dict.lookup("fields") >> selectFields_;
    wordRes::inplaceUniq(selectFields_);

    // Actions to define selection
    selection_ = dict.subOrEmptyDict("selection");

    return true;
}


bool Foam::functionObjects::ensightWrite::execute()
{
    return true;
}

bool Foam::functionObjects::ensightWrite::createCase()
{
    bool staticRestart = false;
    if (ensCase_.valid())
    {
        return staticRestart;
    }

    const Time& t = obr_.time();

    // Define sub-directory name to use for EnSight data.
    // The path to the ensight directory is at case level only
    // - For parallel cases, data only written from master

    fileName ensightDir = dirName_/name();
    if (!ensightDir.isAbsolute())
    {
        ensightDir = t.rootPath()/t.globalCaseName()/ensightDir;
    }

    //Check if case file already exists
    fileName caseFile(ensightDir/t.globalCaseName()+".case");

    DynamicList<scalar> timeSteps(100);
    DynamicList<label> timeIndexes(100);
    if (isFile(caseFile))
    {
        IFstream is(caseFile);
        bool timestepActive = false;
        bool indexActive = false;
        label nTimeStep = 0;
        label nSet = 0;

        //parse existing case file and find output times
        while (is.good())
        {
            string line;

            is.getLine(line);

            // break conditions
            if (line.substr(0, 5) == "# end")
            {
                break;
            }
            else if (line.substr(0, 9) == "time set:")
            {
                IStringStream lineStream(line);
                word dummy;
                label timeSet;
                lineStream >> dummy >> dummy >> timeSet;

                if (timeSet > 1)
                {
                    break;
                }
            }

            if (line.substr(0, 16) == "number of steps:")
            {
                IStringStream lineStream(line);
                word dummy;
                lineStream >> dummy >> dummy >> dummy >> nTimeStep;
            }

            if (indexActive)
            {
                IStringStream lineStream(line);
                for (label j=0; j < 6; j++)
                {
                    label timeIndex;
                    if (nSet < nTimeStep)
                    {
                        lineStream >> timeIndex;
                        timeIndexes.append(timeIndex);
                        nSet++;
                    }
                }
            }

            if (timestepActive)
            {
                IStringStream lineStream(line);
                for (label j=0; j < 6; j++)
                {
                    scalar timeVal;
                    if (nSet < nTimeStep)
                    {
                        lineStream >> timeVal;
                        timeSteps.append(timeVal);
                        nSet++;
                    }
                }
            }

            if (line.substr(0, 12) == "time values:")
            {
                nSet = 0;
                timestepActive = true;
                indexActive = false;
            }

            if (line.substr(0, 17) == "filename numbers:")
            {
                nSet = 0;
                indexActive = true;
            }
        }
    }

    // Wait for all threads to get here.
    Pstream::barrier();

    ensCase_.reset
    (
        new ensightCase
        (
            ensightDir,
            t.globalCaseName(),
            caseOpts_
         )
    );

    //Reset field and geometry output times after a restart
    if (timeSteps.size() > 0)
    {
        DynamicList<label> geomTimes(timeSteps.size());
        forAll(timeSteps, i)
        {
            scalar timeVal = timeSteps[i];

            if (t.value() > timeVal)
            {
                label timeIndex = -1;
                fileName geomFile;
                if (consecutive_)
                {
                    timeIndex = i;
                    geomFile =
                        ensightDir/"data"/ensCase().padded(timeIndex)/"geometry";
                    ensCase().nextTime(timeVal);
                }
                else
                {
                    timeIndex = timeIndexes[i];
                    geomFile =
                        ensightDir/"data"/ensCase().padded(timeIndex)/"geometry";
                    ensCase().setTime(timeVal, timeIndex);
                }
                if (isFile(geomFile))
                {
                    geomTimes.append(timeIndex);
                }
            }
        }

        //Perform sync to ensure all threads are read before master removal
        Pstream::barrier();

        forAll(timeSteps, i)
        {
            scalar timeVal = timeSteps[i];
            if (t.value() <= timeVal)
            {
                label timeIndex = -1;
                if (consecutive_)
                {
                    timeIndex = i;
                }
                else
                {
                    timeIndex = timeIndexes[i];
                }
                fileName timeDir = ensightDir/"data"/ensCase().padded(timeIndex);
                if (isDir(timeDir) && Pstream::master())
                {
                    rmDir(timeDir);
                }
            }
        }

        if (geomTimes.size() > 0)
        {
            if (geomTimes.size() == 1 && geomTimes.size() < timeSteps.size())
            {
                staticRestart = true;
            }

            ensCase().setGeomTimes(geomTimes);
        }
    }
    return staticRestart;
}


bool Foam::functionObjects::ensightWrite::write()
{
    const Time& t = obr_.time();

    bool staticRestart = createCase();

    if (consecutive_)
    {
        ensCase().nextTime(t.value());
    }
    else
    {
        ensCase().setTime(t.value(), t.timeIndex());
    }

    if (update())
    {
        if (!staticRestart)
        {
            autoPtr<ensightGeoFile> os = ensCase_().newGeometry(true);
            ensMesh_().write(os);
        }
    }

    Log << type() << " " << name() << " write: (";

    wordHashSet candidates(subsetStrings(selectFields_, mesh_.names()));
    DynamicList<word> missing(selectFields_.size());
    DynamicList<word> ignored(selectFields_.size());

    // check exact matches first
    forAll(selectFields_, i)
    {
        const wordRe& select = selectFields_[i];
        if (!select.isPattern())
        {
            const word& fieldName = static_cast<const word&>(select);

            if (!candidates.erase(fieldName))
            {
                missing.append(fieldName);
            }
            else if (process(fieldName) < 1)
            {
                ignored.append(fieldName);
            }
        }
    }

    forAllConstIter(wordHashSet, candidates, iter)
    {
        process(iter.key());
    }

    Log << " )" << endl;

    if (missing.size())
    {
        WarningInFunction
            << "Missing field " << missing << endl;
    }
    if (ignored.size())
    {
        WarningInFunction
            << "Unprocessed field " << ignored << endl;
    }

    //Write the case file
    ensCase().write();

    return true;
}


bool Foam::functionObjects::ensightWrite::end()
{
    if (ensCase_.valid())
    {
        // finalize case
        ensCase().write();
        ensCase_.clear();
    }

    return true;
}

// ************************************************************************* //
