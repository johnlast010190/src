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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/functionObjects/writeFile/writeFile.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionObjects::writeFile::outputPrefix
(
    "postProcessing"
);

Foam::label Foam::functionObjects::writeFile::addChars = 7;


const Foam::Enum
<
    Foam::functionObjects::writeFile::formats
>
Foam::functionObjects::writeFile::formatNames_
{
    {Foam::functionObjects::writeFile::formats::dat, "dat"},
    {Foam::functionObjects::writeFile::formats::csv, "csv"},
    {Foam::functionObjects::writeFile::formats::tsv, "tsv"}
};



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeFile::initStream(Ostream& os) const
{
    os.setf(ios_base::scientific, ios_base::floatfield);
    os.precision(writePrecision_);
    os.width(charWidth());
}


Foam::fileName Foam::functionObjects::writeFile::baseFileDir() const
{
    fileName baseDir = fileObr_.time().path();

    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        baseDir = baseDir/".."/outputPrefix;
    }
    else
    {
        baseDir = baseDir/outputPrefix;
    }

    // Append mesh name if not default region
    if (isA<polyMesh>(fileObr_))
    {
        const polyMesh& mesh = refCast<const polyMesh>(fileObr_);
        if (mesh.name() != polyMesh::defaultRegion)
        {
            baseDir = baseDir/mesh.name();
        }
    }
    else if (isA<polyMesh>(fileObr_.parent()))
    {
        // Append mesh and region name if different
        const polyMesh& mesh = refCast<const polyMesh>(fileObr_.parent());
        if (mesh.name() != polyMesh::defaultRegion)
        {
            baseDir = baseDir/mesh.name();
        }
        if (fileObr_.name() != word::null && fileObr_.name() != mesh.name())
        {
            baseDir = baseDir/fileObr_.name();
        }
    }

    // Remove any ".."
    baseDir.clean();

    return baseDir;
}


Foam::fileName Foam::functionObjects::writeFile::baseTimeDir() const
{
    return baseFileDir()/prefix_/fileObr_.time().timeName();
}


void Foam::functionObjects::writeFile::shiftPrevFile
(
    const fileName& outputDir,
    const word fName,
    label& level
) const
{
    label maxLevel = 10;

    word localFname = fName;

    if (level != 0)
    {
        localFname = fName + "_" + name(level - 1);
    }

    fileName sourceFile(outputDir/(localFname + "." + formatNames_[format_]));

    IFstream is(sourceFile);
    if (is.good() && level < maxLevel)
    {
        level++;
        shiftPrevFile(outputDir, fName, level);
        level--;
    }

    fileName targetFile(outputDir/(fName + "_" + name(level) + "." + formatNames_[format_]));
    mv(sourceFile, targetFile);

}


Foam::autoPtr<Foam::OFstream> Foam::functionObjects::writeFile::createFile
(
    const word& name
) const
{
    autoPtr<OFstream> osPtr;

    if (Pstream::master() && writeToFile_)
    {
        const scalar startTime = fileObr_.time().startTime().value();
        const scalar userStartTime = fileObr_.time().timeToUserTime(startTime);
        const word startTimeName = Time::timeName(userStartTime);

        fileName outputDir(baseFileDir()/prefix_/startTimeName);

        mkDir(outputDir);

        word fName(name);

        // Shift existing files to make place for new output file
        label level = 0;
        shiftPrevFile(outputDir, fName, level);

        fileName filePath = outputDir/(fName + "." + formatNames_[format_]);

#if defined(WIN64) || defined(WIN32)
        if (filePath.size() > 256)
        {
            WarningInFunction << "The following function object file path "
                << " exceeds maximum allowable 256 characters for MS Windows : "
                << filePath << nl
                << endl;
        }
#endif

        osPtr.set(new OFstream(filePath));

        initStream(osPtr());
    }

    return osPtr;
}


void Foam::functionObjects::writeFile::resetFile(const word& fileName)
{
    fileName_ = fileName;
    filePtr_ = createFile(fileName_);
}


Foam::Omanip<int> Foam::functionObjects::writeFile::valueWidth
(
    const label offset
) const
{
    return setw(writePrecision_ + addChars + offset);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const word& prefix
)
:
    fileObr_(obr),
    prefix_(prefix),
    fileName_("undefined"),
    filePtr_(),
    writePrecision_(IOstream::defaultPrecision()),
    writeToFile_(true)
{}


Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const word& prefix,
    const word& fileName,
    const dictionary& dict
)
:
    fileObr_(obr),
    prefix_(prefix),
    fileName_(fileName),
    filePtr_(),
    writePrecision_(IOstream::defaultPrecision()),
    writeToFile_(true)
{
    read(dict);

    if (writeToFile_)
    {
        filePtr_ = createFile(fileName_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeFile::~writeFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeFile::read(const dictionary& dict)
{
    writePrecision_ =
        dict.lookupOrDefault("writePrecision", IOstream::defaultPrecision());

    // Only write on master process
    writeToFile_ = dict.lookupOrDefault("writeToFile", true);
    writeToFile_ = writeToFile_ && Pstream::master();

    format_ = formatNames_.lookupOrDefault("outputFileFormat", dict, formats::dat);
    switch(format_)
    {
        case formats::dat:
            delimiter_ = token::TAB;
            spacesAfterDelimiter_ = 0;
            setWidth_ = true;
            break;
        case formats::csv:
            delimiter_ = token::COMMA;
            spacesAfterDelimiter_ = 1;
            setWidth_ = false;
            break;
        case formats::tsv:
            delimiter_ = token::TAB;
            spacesAfterDelimiter_ = 0;
            setWidth_ = false;
            break;
        default:
            delimiter_ = token::TAB;
            spacesAfterDelimiter_ = 0;
            setWidth_ = true;
    }

    return true;
}


Foam::OFstream& Foam::functionObjects::writeFile::file()
{
    if (!writeToFile_)
    {
        return Snull;
    }

    if (!filePtr_.valid())
    {
        FatalErrorInFunction
            << "File pointer not allocated";
    }

    return filePtr_();
}


bool Foam::functionObjects::writeFile::writeToFile() const
{
    return writeToFile_;
}


Foam::label Foam::functionObjects::writeFile::charWidth() const
{
    return writePrecision_ + addChars;
}


void Foam::functionObjects::writeFile::writeCommented
(
    Ostream& os,
    const string& str
) const
{
    unsigned short width = 0;
    if (setWidth_)
    {
        width = charWidth() - 2;
    }
    os  << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(width) << str.c_str();
}

void Foam::functionObjects::writeFile::writeTabbed
(
    Ostream& os,
    const string& str
) const
{
    writeDelimited(os, str);
}

void Foam::functionObjects::writeFile::writeHeader
(
    Ostream& os,
    const string& str
) const
{
    writeCommented(os, str);
    os  << nl;
}


void Foam::functionObjects::writeFile::writeTime(Ostream& os) const
{
    unsigned short width = 0;
    if (setWidth_)
    {
        width = charWidth();
    }
    const scalar timeNow = fileObr_.time().timeOutputValue();
    os  << setw(width) << Time::timeName(timeNow);
}


const Foam::word Foam::functionObjects::writeFile::getFormattedTime() const
{
    return Time::timeName(fileObr_.time().timeOutputValue());
}


// ************************************************************************* //
