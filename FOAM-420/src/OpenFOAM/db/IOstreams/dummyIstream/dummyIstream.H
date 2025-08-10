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
    (c) 2017 OpenFOAM Foundation

Class
    Foam::dummyIstream

Description
    Dummy stream for input. Aborts at any attempt to read from it.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef dummyIstream_H
#define dummyIstream_H

#include "db/IOstreams/IOstreams/Istream.H"
//#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dummyIstream Declaration
\*---------------------------------------------------------------------------*/

class dummyIstream
:
    public Istream
{

public:

    // Constructors

        //- Construct null
        dummyIstream()
        :
            Istream()
        {}


    //- Destructor
    ~dummyIstream()
    {}


    // Member functions

        // Read functions

            //- Return next token from stream
            virtual Istream& read(token&)
            {
                NotImplemented;
            }

            //- Read a character
            virtual Istream& read(char&)
            {
                NotImplemented;
            }

            //- Read a word
            virtual Istream& read(word&)
            {
                NotImplemented;
            }

            // Read a string (including enclosing double-quotes)
            virtual Istream& read(string&)
            {
                NotImplemented;
            }

            //- Read a label
            virtual Istream& read(label&)
            {
                NotImplemented;
            }

            //- Read a floatScalar
            virtual Istream& read(floatScalar&)
            {
                NotImplemented;
            }

            //- Read a doubleScalar
            virtual Istream& read(doubleScalar&)
            {
                NotImplemented;
            }

            //- Read binary block
            virtual Istream& read(char*, std::streamsize)
            {
                NotImplemented;
            }

            //- Rewind and return the stream so that it may be read again
            virtual Istream& rewind()
            {
                NotImplemented;
            }

            //- Return flags of stream
            virtual ios_base::fmtflags flags() const
            {
                NotImplemented;
            }

            //- Set flags of stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags f)
            {
                NotImplemented;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
