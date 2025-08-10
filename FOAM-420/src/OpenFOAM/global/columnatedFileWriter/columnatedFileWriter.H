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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::columnatedFileWriter

Description
    Class for writing columnated files (e.g. .dat files, .tsv files, .csv files)

SourceFiles
    columnatedFileWriter.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/objectRegistry/objectRegistry.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "db/IOstreams/token/token.H"
#include "primitives/enums/Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class columnatedFileWriter Declaration
\*---------------------------------------------------------------------------*/

class columnatedFileWriter
{

protected:

    // Protected data

        //- Enum of the possible output types, also used as file extensions
        enum formats
        {
            dat,
            csv,
            tsv
        };

        static const Enum<formats> formatNames_;

        //- Format specified in the input dictionary
        const formats format_ = formats::dat;

        //- Reference to time (used to write out the time to file)
        const Time& time_;

        //- The directory in which to write the output file
        const fileName outputDir_;

        //- The name of the output file (without the path to the file or the
        //  file extension)
        const word outFileName_;

        //- Output file
        OFstream file_;

        //- Ostream for file (for convenience)
        OSstream& os_;

        //- Write precision
        label writePrecision_;

        //- Delimiter
        token delimiter_ = token::TAB;

        //- Number of spaces after the delimiter
        unsigned int spacesAfterDelimiter_ = 0;

        //- Whether to print set-width columns or not
        bool fixedWidth_ = true;


    // Protected Member Functions

        //- Create the output file and return an autoPtr to the OFstream for
        //  writing to the new file.
        //- This is heavily dependent on initialisation order:
        // -    file_ mustn't be initialised yet, because then shiftPrevFiles()
        //      will shift the file we want to write to,
        // -    shiftPrevFiles() required outputDir_, formatNames_, format_, and
        //      outFileName_ to be initialised.
        void prepareOutputDirectory();

        //- Return the value width when writing to stream with optional offset
        Omanip<int> valueWidth(const label offset = 0) const;

        //- Backup old result files by recursively moving them to larger
        //  numbered indices. Only shifts the first contiguous block of existing
        //  files, starting from the original file name.  For example, for a
        //  file U, if files U, U_1, U_3, and U_4 exist, U_1->U_2, U->U_1
        //  (leaving space to write a new file called U), and U_3 and U_4 are
        //  untouched.
        void shiftPrevFiles() const;


public:

    //- Directory prefix
    static const word outputPrefix;

    //- Additional characters for writing
    static label addChars;


    // Constructors

        //- Construct from objectRegistry and prefix, and read options
        //  from dictionary
        columnatedFileWriter
        (
            const Time& time,
            const fileName outputDir,
            const word& outFileName,
            const dictionary& dict
        );


    //- Destructor

        virtual ~columnatedFileWriter() = default;


    // Member Functions

        //- Direct access to the output file stream
        //- Use ONLY if your output will be compatible with .tsv, .csv, and
        //  legacy .dat file formats
        virtual OFstream& file();

        //- True if writing a .dat file.  Can be useful to define specific rules
        //  to maintain legacy file formats.
        bool isDatFile() const;

        //- Return width of character stream output
        virtual label charWidth() const;

        //- Write a commented string to stream.  For column names, prefer
        //  writeDelimitedComment().
        virtual void writeCommented(const string& str) const;

        //- Write a delimited row, but with a "# " preprended to the front.
        //  Useful when writing headers.
        template<typename... Args>
        void writeDelimitedComment
        (
            const string& str,
            const Args&... args
        ) const;

        //- Write a delimited row, but with a "# " preprended to the front.
        //  Useful when writing headers.
        template<typename T1, typename... Args>
        void writeDelimitedComment
        (
            const T1& value,
            const Args&... args
        ) const;

        //- A template used to output columnated data.  This is less useful than
        //  the variadic version of the function, because you have to call this
        //  once for every column of data.
        template<typename T1>
        void writeDelimited(const T1& value) const;

        //- A variadic template allowing the specification of a whole delimited
        //  row all at once.  Use for columnated data.
        template<typename T1, typename... Args>
        void
        writeDelimited(const T1& value, const Args&... args) const;

        //- Write the current time to stream.  Doesn't deal with delimiters, so
        //  assumes Time is the first column.
        virtual void writeTime() const;

        //- Get a nicely-formatted version of the current time as a word.
        const Foam::word getFormattedTime() const;

        //- Write a new line to file
        void newLine() const;

        //- Write a new line to file and flush the output
        void endLine() const;

        //- Write a (commented) header property and value pair
        template<class Type>
        void writeHeaderValue
        (
            const string& property,
            const Type& value
        ) const;

        //- Manually set a delimiter
        //- Useful for legacy formats that don't require special handling, but
        //  don't comply with the setup in the columnatedFileWriter constructor
        void setDelimiter(const token delim);

        //- Manually set the number of spaces after the delimiter
        //- Useful for legacy formats that don't require special handling, but
        //  don't comply with the setup in the columnatedFileWriter constructor
        void setNumberOfSpacesAfterDelimiter(const label spacesAfterDelim);

        //- Manually toggle fixed column width
        //- Useful for legacy formats that don't require special handling, but
        //  don't comply with the setup in the columnatedFileWriter constructor
        void setFixedWidth(bool fixedWidthSwitch);
};


//- Specialisations for strings (don't print "" around strings, for
//  backwards-compatibiliy)
template<>
inline void columnatedFileWriter::writeDelimited<string>(const string& str) const;
template<>
inline void columnatedFileWriter::writeDelimited<std::string>(const std::string& str) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "global/columnatedFileWriter/columnatedFileWriterTemplates.C"
#endif

// ************************************************************************* //
