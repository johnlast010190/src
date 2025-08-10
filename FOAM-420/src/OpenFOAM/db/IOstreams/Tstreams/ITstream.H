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
    (c) 2017-2022 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::ITstream

Description
    Input token stream.

SourceFiles
    ITstream.C

\*---------------------------------------------------------------------------*/

#ifndef ITstream_H
#define ITstream_H

#include "db/IOstreams/IOstreams/Istream.H"
#include "db/IOstreams/token/tokenList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ITstream Declaration
\*---------------------------------------------------------------------------*/

class ITstream
:
    public Istream,
    public tokenList
{
    // Private data

        //- Name of ITstream
        fileName name_;

        //- Index of token currently being read
        label tokenIndex_;

        //- Failsafe read-access to token at specified location
        //- or undefinedToken
        inline const token& peekAt(const label i) const
        {
            return
            (
                i >= 0 && i < tokenList::size()
                ? tokenList::operator[](i)
                : token::undefinedToken
            );
        }


public:

    // Constructors

        //- Construct from components
        ITstream
        (
            const string& name,
            const UList<token>& tokens,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            Istream(format, version),
            tokenList(tokens),
            name_(name),
            tokenIndex_(0)
        {
            setOpened();
            setGood();
        }


        //- Construct from components, transferring the tokens
        ITstream
        (
            const string& name,
            const Xfer<List<token>>& tokens,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            Istream(format, version),
            tokenList(tokens),
            name_(name),
            tokenIndex_(0)
        {
            setOpened();
            setGood();
        }


        //- Construct as copy
        ITstream(const ITstream& its)
        :
            Istream(ASCII, currentVersion),
            tokenList(its),
            name_(its.name_),
            tokenIndex_(0)
        {
            setOpened();
            setGood();
        }


        //- Destructor
        virtual ~ITstream()
        {}


    // Member functions

        // Inquiry

            //- Return the name of the stream
            const fileName& name() const
            {
                return name_;
            }

            //- Return non-const access to the name of the stream
            fileName& name()
            {
                return name_;
            }

    // Token Access

            //- True if putback token is in use
            bool hasPutback() const noexcept
            {
                return Istream::hasPutback();
            }

            //- Failsafe peek at the \b first token in the list.
            //  \return \c undefinedToken if the list is empty.
            const token& front() const { return peekAt(0); }

            //- Failsafe peek at the \b last token in the list.
            //  \return \c undefinedToken if the list is empty.
            const token& back() const { return peekAt(tokenList::size()-1); }

            //- Failsafe peek at what the next read would return,
            //- including handling of any putback
            //  \return \c undefinedToken if list is exhausted
            const token& peek() const;

            //- Return the current token index
            label tokenIndex() const
            {
                return tokenIndex_;
            }

            //- Return non-const access to the current token index
            label& tokenIndex()
            {
                return tokenIndex_;
            }

            //- Return the number of remaining tokens
            label nRemainingTokens() const
            {
                return size() - tokenIndex_;
            }

            //- Return flags of output stream
            ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Read functions

            //- Return next token from stream
            virtual Istream& read(token&);

            //- Read a character
            virtual Istream& read(char&);

            //- Read a word
            virtual Istream& read(word&);

            // Read a string (including enclosing double-quotes)
            virtual Istream& read(string&);

            //- Read a label
            virtual Istream& read(label&);

            //- Read a floatScalar
            virtual Istream& read(floatScalar&);

            //- Read a doubleScalar
            virtual Istream& read(doubleScalar&);

            //- Read binary block
            virtual Istream& read(char*, std::streamsize);

            //- Rewind and return the stream so that it may be read again
            virtual Istream& rewind();


        // Edit

            //- Set flags of stream
            ios_base::fmtflags flags(const ios_base::fmtflags)
            {
                return ios_base::fmtflags(0);
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;


    // Member Operators

        //- Copy assignment, with rewind()
        void operator=(const ITstream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
