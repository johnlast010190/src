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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2022 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"
#include "db/IOstreams/StringStreams/StringStream.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/dictionary/dictionary.H"
#include "global/JobInfo/JobInfo.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "include/OSspecific.H"
#include "algorithms/serialThreads/serialThreads.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::Ostream& Foam::error::printDeprecationWarning
(
    Ostream& os,
    const word& typeName,
    const char* what,
    const int oldVersion
)
{
    if (oldVersion < 10000)
    {
        // 4 digits - something from OPENFOAM+, so merged in FOAM v3
        // Emit warning
        os  << typeName << ": This " << what
            << " was deprecated in FOAM version 3.0.0 and will be removed"
            << " in future.\n"
            << endl << "   ";
    }
    else if (FOAM_API >= oldVersion)
    {
        os  << typeName << ": This " << what
            << " was deprecated in FOAM version "
            << oldVersion/10000 << "."
            << (oldVersion/100) % 100 << "."
            << oldVersion % 100
            << " and will be removed in future."
            << endl << "    ";
    }
    return os;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::error::error(const string& title)
:
    std::exception(),
    messageStream(title, messageStream::FATAL),
    functionName_("unknown"),
    sourceFileName_("unknown"),
    sourceFileLineNumber_(0),
    throwExceptions_(false),
    messageStreamPtr_(new OStringStream())
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::error(const string& title) : cannot open error stream"
            << endl;
        exit(1);
    }
}


Foam::error::error(const dictionary& errDict)
:
    std::exception(),
    messageStream(errDict),
    functionName_(errDict.lookup("functionName")),
    sourceFileName_(errDict.lookup("sourceFileName")),
    sourceFileLineNumber_(readLabel(errDict.lookup("sourceFileLineNumber"))),
    throwExceptions_(false),
    messageStreamPtr_(new OStringStream())
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::error(const dictionary& errDict) : "
               "cannot open error stream"
            << endl;
        exit(1);
    }
}


Foam::error::error(const error& err)
:
    std::exception(),
    messageStream(err),
    functionName_(err.functionName_),
    sourceFileName_(err.sourceFileName_),
    sourceFileLineNumber_(err.sourceFileLineNumber_),
    throwExceptions_(err.throwExceptions_),
    messageStreamPtr_(new OStringStream(*err.messageStreamPtr_))
{
    //*messageStreamPtr_ << err.message();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::error::~error() throw()
{
    delete messageStreamPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::OSstream& Foam::error::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    functionName_ = functionName;
    sourceFileName_ = sourceFileName;
    sourceFileLineNumber_ = sourceFileLineNumber;

    return operator OSstream&();
}


Foam::OSstream& Foam::error::operator()
(
    const string& functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    return operator()
    (
        functionName.c_str(),
        sourceFileName,
        sourceFileLineNumber
    );
}


Foam::error::operator Foam::OSstream&()
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::operator OSstream&() : error stream has failed"
            << endl;
        abort();
    }

    return *messageStreamPtr_;
}


Foam::error::operator Foam::dictionary() const
{
    dictionary errDict;

    string oneLineMessage(message());
    oneLineMessage.replaceAll('\n', ' ');

    errDict.add("type", word("Foam::error"));
    errDict.add("message", oneLineMessage);
    errDict.add("function", functionName());
    errDict.add("sourceFile", sourceFileName());
    errDict.add("sourceFileLineNumber", sourceFileLineNumber());

    return errDict;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::string Foam::error::message() const
{
    return messageStreamPtr_->str();
}


void Foam::error::exit(const int errNo)
{
    if (!throwExceptions_ && JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
        jobInfo.exit();
    }

    if (env("FOAM_ABORT") || env("FOAM_ABORT"))
    {
        abort();
    }

    if (throwExceptions_)
    {
        // Make a copy of the error to throw
        error errorException(*this);

        // Reset the message buffer for the next error message
        messageStreamPtr_->reset();

        throw errorException;
    }
    else if (Pstream::parRun())
    {
        Perr<< endl << *this << endl
            << "\nFOAM parallel run exiting\n" << endl;
        serialThreads::end();
        Pstream::exit(errNo);
    }
    else
    {
        Perr<< endl << *this << endl
            << "\nFOAM exiting\n" << endl;
        serialThreads::end();
        ::exit(1);
    }
}


void Foam::error::abort()
{
    if (!throwExceptions_ && JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
        jobInfo.abort();
    }

    if (env("FOAM_ABORT") || env("FOAM_ABORT"))
    {
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT or FOAM_ABORT set)\n" << endl;
        printStack(Perr);
        serialThreads::end();
        ::abort();
    }

    if (throwExceptions_)
    {
        // Make a copy of the error to throw
        error errorException(*this);

        // Reset the message buffer for the next error message
        messageStreamPtr_->reset();

        throw errorException;
    }

    if (Pstream::parRun())
    {
        Perr<< endl << *this << endl
            << "\nFOAM parallel run aborting\n" << endl;
        printStack(Perr);
        serialThreads::end();
        Pstream::abort();
    }
    else
    {
        Perr<< endl << *this << endl
            << "\nFOAM aborting\n" << endl;
        printStack(Perr);
        serialThreads::end();
        ::abort();
    }
}


void Foam::error::write(Ostream& os, const bool includeTitle) const
{
    os  << nl;
    if (includeTitle)
    {
        os  << title().c_str() << endl;
    }
    os  << message().c_str();

    if (error::level >= 2 && sourceFileLineNumber())
    {
        os  << nl << nl
            << "    From function " << functionName().c_str() << endl
            << "    in file " << sourceFileName().c_str()
            << " at line " << sourceFileLineNumber() << '.';
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const error& err)
{
    err.write(os);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

Foam::error Foam::FatalError("--> FOAM FATAL ERROR: ");


// ************************************************************************* //
