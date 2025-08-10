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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2011 Symscape
    (c) 2014 blueCAPE Lda

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.



\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Fstreams/OFstream.H"
#include "include/OSspecific.H"
#include "db/IOstreams/gzstream/gzstream.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstream, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::OFstreamAllocator::OFstreamAllocator
(
    const fileName& pathname,
    IOstream::compressionType compression,
    const bool append
)
:
    allocatedPtr_(nullptr)
{
    if (pathname.empty())
    {
        if (OFstream::debug)
        {
            InfoInFunction << "Cannot open null file " << endl;
        }
    }
    ofstream::openmode mode(ofstream::out);
    if (append)
    {
        mode |= ofstream::app;
    }

    if (compression == IOstream::COMPRESSED)
    {
        // Get identically named uncompressed version out of the way
        fileName::Type pathType = Foam::type(pathname, false);
        if (pathType == fileName::FILE || pathType == fileName::LINK)
        {
            rm(pathname);
        }
        fileName gzPathName(pathname + ".gz");

        if (!append && Foam::type(gzPathName) == fileName::LINK)
        {
            // Disallow writing into softlink to avoid any problems with
            // e.g. softlinked initial fields
            rm(gzPathName);
        }

        allocatedPtr_ = new ogzstream(gzPathName.c_str(), mode);
    }
    else
    {
        // get identically named compressed version out of the way
        fileName gzPathName(pathname + ".gz");
        fileName::Type gzType = Foam::type(gzPathName, false);
        if (gzType == fileName::FILE || gzType == fileName::LINK)
        {
            rm(gzPathName);
        }
        if (!append && Foam::type(pathname, false) == fileName::LINK)
        {
            // Disallow writing into softlink to avoid any problems with
            // e.g. softlinked initial fields
            rm(pathname);
        }

#if defined(WIN64) || defined(WIN32)
        // Use binary mode in case we write binary.
        // Causes windows reading to fail if we don't
        allocatedPtr_ = new ofstream(pathname.c_str(),
                              ios_base::out|ios_base::binary);
#else
        allocatedPtr_ = new ofstream(pathname.c_str(),mode);
#endif
    }
}


Foam::OFstreamAllocator::~OFstreamAllocator()
{
    deallocate();
}


void Foam::OFstreamAllocator::deallocate()
{
    if (allocatedPtr_)
    {
        delete allocatedPtr_;
        allocatedPtr_ = nullptr;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstream::OFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version,
    compressionType compression,
    const bool append
)
:
    OFstreamAllocator(pathname, compression, append),
    OSstream
    (
        *allocatedPtr_,
        "OFstream.sinkFile_",
        format,
        version,
        compression
    )
{
    setClosed();
    setState(allocatedPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for output" << nl << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstream::~OFstream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::ostream& Foam::OFstream::stdStream()
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *allocatedPtr_;
}


const std::ostream& Foam::OFstream::stdStream() const
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *allocatedPtr_;
}


void Foam::OFstream::print(Ostream& os) const
{
    os  << "OFstream: ";
    OSstream::print(os);
}


// ************************************************************************* //
