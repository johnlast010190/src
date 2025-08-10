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
    (c) 2011 Symscape
    (c) 2014 blueCAPE Lda

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.



\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Fstreams/IFstream.H"
#include "include/OSspecific.H"
#include "db/IOstreams/gzstream/gzstream.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IFstream, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::IFstreamAllocator::IFstreamAllocator(const fileName& pathname)
:
    ifPtr_(nullptr),
    compression_(IOstream::UNCOMPRESSED)
{
    if (pathname.empty())
    {
        if (IFstream::debug)
        {
            InfoInFunction << "Cannot open null file " << endl;
        }
    }

#if defined(WIN64) || defined(WIN32)
    // Use binary mode in case we read binary.
    // Causes windows reading to fail if we don't.
    ifPtr_ = new ifstream(pathname.c_str(),
                          ios_base::in|ios_base::binary);
#else
    ifPtr_ = new ifstream(pathname.c_str());
#endif

    // If the file is compressed, decompress it before reading.
    if (!ifPtr_->good() && isFile(pathname + ".gz", false))
    {
        if (IFstream::debug)
        {
            InfoInFunction << "Decompressing " << pathname + ".gz" << endl;
        }

        delete ifPtr_;

        ifPtr_ = new igzstream((pathname + ".gz").c_str());

        if (ifPtr_->good())
        {
            compression_ = IOstream::COMPRESSED;
        }
    }
}


Foam::IFstreamAllocator::~IFstreamAllocator()
{
    delete ifPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IFstream::IFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version
)
:
    IFstreamAllocator(pathname),
    ISstream
    (
        *ifPtr_,
        "IFstream.sourceFile_",
        format,
        version,
        IFstreamAllocator::compression_
    ),
    pathname_(pathname)
{
    setClosed();

    setState(ifPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            InfoInFunction
                << "Could not open file for input" << endl << info() << endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::istream& Foam::IFstream::stdStream()
{
    if (!ifPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated" << abort(FatalError);
    }
    return *ifPtr_;
}


const std::istream& Foam::IFstream::stdStream() const
{
    if (!ifPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated" << abort(FatalError);
    }
    return *ifPtr_;
}


void Foam::IFstream::print(Ostream& os) const
{
    // Print File data
    os  << "IFstream: ";
    ISstream::print(os);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::IFstream& Foam::IFstream::operator()() const
{
    if (!good())
    {
        // also checks .gz file
        if (isFile(pathname_, true))
        {
            check("IFstream::operator()");
            FatalIOError.exit();
        }
        else
        {
            FatalIOErrorInFunction(*this)
                << "file " << pathname_ << " does not exist"
                << exit(FatalIOError);
        }
    }

    return const_cast<IFstream&>(*this);
}


// ************************************************************************* //
