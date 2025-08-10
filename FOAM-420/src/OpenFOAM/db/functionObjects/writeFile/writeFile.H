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
    Foam::functionObjects::writeFile

Description
    functionObject base class for writing single files

See also
    Foam::functionObject
    Foam::functionObjects::logFiles

SourceFiles
    writeFile.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeFile_H
#define functionObjects_writeFile_H

#include "db/objectRegistry/objectRegistry.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "db/IOstreams/token/token.H"
#include "primitives/enums/Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class writeFile Declaration
\*---------------------------------------------------------------------------*/

class writeFile
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

        //- Format specified on command-line
        formats format_ = formats::dat;

        //- Reference to the region objectRegistry
        const objectRegistry& fileObr_;

        //- Prefix
        const word prefix_;

        //- Name of file
        word fileName_;

        //- File pointer
        autoPtr<OFstream> filePtr_;

        //- Write precision
        label writePrecision_;

        //- Flag to enable/disable writing to file
        bool writeToFile_;


        //- Delimiter
        token delimiter_ = token::TAB;

        //- Number of spaces after the delimiter
        unsigned short spacesAfterDelimiter_ = 0;

        //- Whether to print set-width columns or not
        bool setWidth_ = true;


    // Protected Member Functions

        //- Initialise the output stream for writing
        void initStream(Ostream& os) const;

        //- Return an autoPtr to a new file
        virtual autoPtr<OFstream> createFile(const word& name) const;

        //- Reset internal file pointer to new file with new name
        virtual void resetFile(const word& name);

        //- Return the value width when writing to stream with optional offset
        Omanip<int> valueWidth(const label offset = 0) const;

        //- backup old result files by recursively moving them to larger
        //numbered indices
        void shiftPrevFile(const fileName&, const word, label&) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        writeFile(const writeFile&);

        //- Disallow default bitwise assignment
        void operator=(const writeFile&);


public:

    //- Directory prefix
    static const word outputPrefix;

    //- Additional characters for writing
    static label addChars;


    // Constructors

        //- Construct from objectRegistry and prefix
        writeFile
        (
            const objectRegistry& obr,
            const word& prefix
        );

        //- Construct from objectRegistry and prefix, and read options
        //  from dictionary
        writeFile
        (
            const objectRegistry& obr,
            const word& prefix,
            const word& fileName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeFile();


    // Member Functions

        //- Read
        virtual bool read(const dictionary& dict);

        //- Return access to the file (if only 1)
        virtual OFstream& file();

        //- Flag to allow writing to file
        virtual bool writeToFile() const;

        //- Return width of character stream output
        virtual label charWidth() const;

        //- Write a commented string to stream.  For column names, prefer
        //  writeDelimitedComment().
        virtual void writeCommented(Ostream& os, const string& str) const;

        //- Write a delimited row, but with a "# " preprended to the front.
        //  Useful when writing headers.
        template<typename... Args>
        void writeDelimitedComment
        (
            Ostream& os,
            const string& str,
            const Args&... args
        ) const;

        //- Write a delimited row, but with a "# " preprended to the front.
        //  Useful when writing headers.
        template<typename T1, typename... Args>
        void writeDelimitedComment
        (
            Ostream& os,
            const T1& value,
            const Args&... args
        ) const;

        //- A template used to output columnated data.  This is less useful than
        //  the variadic version of the function, because you have to call this
        //  once for every column of data.
        template<typename T1>
        void writeDelimited(Ostream& os, const T1& value) const;

        //- A variadic template allowing the specification of a whole delimited
        //  row all at once.  Use for columnated data.
        template<typename T1, typename... Args>
        void
        writeDelimited(Ostream& os, const T1& value, const Args&... args) const;

        //- Write a tabbed string to stream.  writeDelimited() is identical to
        //  writeTabbed if the following conditions are met:
        // - delimiter_ == token::TAB;
        // - spacesAfterDelimiter_ == 0;
        // - setWidth_ == True;
        // * All these conditions are true by default with the "legacy" format.
        // * DEPRECATED:  Prefer writeDelimited().
        FOAM_DEPRECATED_FOR(v3.5, "writeDelimited()")
        virtual void writeTabbed(Ostream& os, const string& str) const;

        //- Write a commented header to stream, i.e. write the provided string,
        //  prepended by the comment symbol, followed by a new line.
        virtual void writeHeader(Ostream& os, const string& str) const;

        //- Write the current time to stream.  Doesn't deal with delimiters, so
        //  assumes Time is the first column.
        virtual void writeTime(Ostream& os) const;

        //- Get a nicely-formatted version of the current time as a word.
        const Foam::word getFormattedTime() const;

        //- Write a (commented) header property and value pair
        template<class Type>
        void writeHeaderValue
        (
            Ostream& os,
            const string& property,
            const Type& value
        ) const;

        //- Return the base directory for output
        fileName baseFileDir() const;

        //- Return the base directory for the current time value
        fileName baseTimeDir() const;
};


//- Specialisations for strings (don't print "" around strings, for
//  backwards-compatibiliy)
template<>
inline void writeFile::writeDelimited<string>(Ostream& os, const string& str) const;
template<>
inline void writeFile::writeDelimited<std::string>(Ostream& os, const std::string& str) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/functionObjects/writeFile/writeFileTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
