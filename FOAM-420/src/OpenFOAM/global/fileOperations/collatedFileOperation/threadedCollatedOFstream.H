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
    (c) 2017-2018 OpenFOAM Foundation

Class
    Foam::threadedCollatedOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    threadedCollatedOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef threadedCollatedOFstream_H
#define threadedCollatedOFstream_H

#include "db/IOstreams/StringStreams/StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstreamCollator;

/*---------------------------------------------------------------------------*\
                  Class threadedCollatedOFstream Declaration
\*---------------------------------------------------------------------------*/

class threadedCollatedOFstream
:
    public OStringStream
{
    // Private data

        OFstreamCollator& writer_;

        const fileName pathName_;

        const IOstream::compressionType compression_;

        const bool useThread_;


public:

    // Constructors

        //- Construct and set stream status
        threadedCollatedOFstream
        (
            OFstreamCollator&,
            const fileName& pathname,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED,
            const bool useThread = true
        );


    //- Destructor
    ~threadedCollatedOFstream();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
