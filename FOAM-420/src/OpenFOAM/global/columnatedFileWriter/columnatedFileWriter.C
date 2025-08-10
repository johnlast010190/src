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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "global/columnatedFileWriter/columnatedFileWriter.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/Fstreams/IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
// TODO:  Review these, looks like some carry-over from old code

const Foam::word Foam::columnatedFileWriter::outputPrefix("postProcessing");

Foam::label Foam::columnatedFileWriter::addChars = 7;


const Foam::Enum<Foam::columnatedFileWriter::formats>
Foam::columnatedFileWriter::formatNames_
{
    {Foam::columnatedFileWriter::formats::dat, "dat"},
    {Foam::columnatedFileWriter::formats::csv, "csv"},
    {Foam::columnatedFileWriter::formats::tsv, "tsv"}
};



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::columnatedFileWriter::shiftPrevFiles() const
{
    auto fullPathToFile = [&] (const word& fName) -> fileName
    {
        return fileName
        (
            outputDir_ /
            (fName + "." + formatNames_[format_])
        );
    };

    // Probably unnecessary, but makes the code nice and readable
    auto fileExists = [&] (const fileName& f) -> bool
    {
        return IFstream(f).good();
    };


    if (!fileExists(fullPathToFile(outFileName_)))
    {
        // There's nothing at the path we intend to write to - early exit!
        return;
    }

    const label maxLevel = 10;
    DynamicList<fileName> filesToMove(maxLevel);
    filesToMove.append(fullPathToFile(outFileName_));

    for (auto i = 1; i < maxLevel; ++i)
    {
        const fileName fileToMove = fullPathToFile(outFileName_ + "_" + name(i));
        if (fileExists(fileToMove))
        {
            filesToMove.append(fileToMove);
        }
        else
        {
            // Can safely move the file at i-1 without overwriting anything
            // (existing files are at least not contiguous)
            break;
        }
    }
    filesToMove.shrink();

    // DynamicList reverse iterator is broken, see CO-6540
    for (auto i = filesToMove.size() - 1; i >= 0; --i)
    {
        mv(filesToMove[i], fullPathToFile(outFileName_ + "_" + name(i + 1)));
    }

}


void Foam::columnatedFileWriter::prepareOutputDirectory()
{
    if (Pstream::master())
    {
        mkDir(outputDir_);
        shiftPrevFiles();
    }
    Pstream::barrier();
}


Foam::Omanip<int>
Foam::columnatedFileWriter::valueWidth(const label offset) const
{
    return setw(writePrecision_ + addChars + offset);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::columnatedFileWriter::columnatedFileWriter
(
    const Time& time,
    const fileName outputDir,
    const word& outFileName,
    const dictionary& dict
)
:
    format_(formatNames_.lookupOrDefault("outputFileFormat", dict, formats::dat)),
    time_(time),
    outputDir_(outputDir),
    outFileName_(outFileName),
    file_
    (
        (
            prepareOutputDirectory(),
            outputDir_/(outFileName_ + "." + formatNames_[format_])
        )
    ),
    os_(Pstream::master() ? file_ : Snull),
    writePrecision_(dict.lookupOrDefault("writePrecision", IOstream::defaultPrecision()))
{
    file_.setf(ios_base::scientific, ios_base::floatfield);
    file_.precision(writePrecision_);
    file_.width(charWidth());

    // TODO:  delimiter_, spacesAfterDelimiter_, and fixedWidth_ should all
    // probably be const, but I can't think of an elegant way to get the
    // construction initialisation working for that
    switch(format_)
    {
        case formats::dat:
            delimiter_ = token::TAB;
            spacesAfterDelimiter_ = 0;
            fixedWidth_ = true;
            break;
        case formats::csv:
            delimiter_ = token::COMMA;
            spacesAfterDelimiter_ = 1;
            fixedWidth_ = false;
            break;
        case formats::tsv:
            delimiter_ = token::TAB;
            spacesAfterDelimiter_ = 0;
            fixedWidth_ = false;
            break;
        default:
            delimiter_ = token::TAB;
            spacesAfterDelimiter_ = 0;
            fixedWidth_ = true;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::OFstream& Foam::columnatedFileWriter::file()
{
    if (not Pstream::master())
    {
        return Snull;
    }

    return file_;
}

bool Foam::columnatedFileWriter::isDatFile() const
{
    return formats::dat == format_;
}


Foam::label Foam::columnatedFileWriter::charWidth() const
{
    return writePrecision_ + addChars;
}


void Foam::columnatedFileWriter::writeCommented
(
    const string& str
) const
{
    unsigned short width = 0;
    if (fixedWidth_)
    {
        width = charWidth() - 2;
    }
    os_ << setw(1) << "#" << setw(1) << ' '
        << setf(ios_base::left) << setw(width) << str.c_str();
}


void Foam::columnatedFileWriter::writeTime() const
{
    unsigned short width = 0;
    if (fixedWidth_)
    {
        width = charWidth();
    }
    const scalar timeNow = time_.timeOutputValue();
    os_ << setw(width) << Time::timeName(timeNow);
}


const Foam::word Foam::columnatedFileWriter::getFormattedTime() const
{
    return Time::timeName(time_.timeOutputValue());
}


void Foam::columnatedFileWriter::newLine() const
{
    os_ << nl;
}


void Foam::columnatedFileWriter::endLine() const
{
    os_ << endl;
}

void Foam::columnatedFileWriter::setDelimiter(const Foam::token delim)
{
    delimiter_ = delim;
}

void Foam::columnatedFileWriter::setNumberOfSpacesAfterDelimiter(const Foam::label spacesAfterDelim)
{
    spacesAfterDelimiter_ = spacesAfterDelim;
}

void Foam::columnatedFileWriter::setFixedWidth(bool fixedWidthSwitch)
{
    fixedWidth_ = fixedWidthSwitch;
}


// ************************************************************************* //
