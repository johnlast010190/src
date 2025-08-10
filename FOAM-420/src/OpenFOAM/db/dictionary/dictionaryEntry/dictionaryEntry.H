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
    (c) 2010-2012 Esi Ltd.
    (c) 2018 OpenCFD Ltd.

Class
    Foam::dictionaryEntry

Description
    A keyword and a list of tokens is a 'dictionaryEntry'.

    An dictionaryEntry can be read, written and printed, and the types and
    values of its tokens analysed.  A dictionaryEntry is a high-level building
    block for data description.  It is a front-end for the token parser.
    A list of entries can be used as a set of keyword syntax elements,
    for example.

SourceFiles
    dictionaryEntry.C
    dictionaryEntryIO.C

\*---------------------------------------------------------------------------*/

#ifndef dictionaryEntry_H
#define dictionaryEntry_H

#include "db/dictionary/entry/entry.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/IOstreams/InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class dictionaryEntry;

Ostream& operator<<(Ostream&, const dictionaryEntry&);


/*---------------------------------------------------------------------------*\
                           Class dictionaryEntry Declaration
\*---------------------------------------------------------------------------*/

class dictionaryEntry
:
    public entry,
    public dictionary
{
    // Private Member Functions

        void readData(Istream&);

        //- Dissallow bitwise copy
       dictionaryEntry(const dictionary&);


public:

    // Constructors

        //- Construct from the parent dictionary and Istream
        dictionaryEntry(const dictionary& parentDict, Istream&);

        //- Construct from the keyword, parent dictionary and a Istream
        dictionaryEntry
        (
            const keyType&,
            const dictionary& parentDict,
            Istream&
        );

        //- Construct from the keyword, parent dictionary and a dictionary
        dictionaryEntry
        (
            const keyType&,
            const dictionary& parentDict,
            const dictionary& dict
        );

        //- Construct as copy for the given parent dictionary
        dictionaryEntry
        (
            const dictionary& parentDict,
            const dictionaryEntry&
        );

        //- construct null
        dictionaryEntry();

        //- Construct as copy
        //dictionaryEntry(const dictionaryEntry&);

        autoPtr<entry> clone(const dictionary& parentDict) const
        {
            return autoPtr<entry>(new dictionaryEntry(parentDict, *this));
        }


    // Member functions

        //- Return the local dictionary name (final part of scoped name)
        const word dictName() const
        {
            return dictionary::dictName();
        }

        //- Return the dictionary name (scoped, e.g. dictA::dictB::dictC)
        const fileName& name() const
        {
            return dictionary::name();
        }

        //- Return the dictionary name (scoped, e.g. dictA::dictB::dictC)
        fileName& name()
        {
            return dictionary::name();
        }

        //- Return line number of first token in dictionary
        label startLineNumber() const;

        //- Return line number of last token in dictionary
        label endLineNumber() const;

        //- This entry is not a primitive,
        //  calling this function generates a FatalError
        ITstream& stream() const;


        //- Return pointer to this dictionary
        virtual const dictionary* dictPtr() const;

        //- Return non-const pointer to this dictionary
        virtual dictionary* dictPtr();

        //- Return true because this entry is a dictionary
        bool isDict() const
        {
            return true;
        }

        //- Return dictionary
        const dictionary& dict() const;

        //- Return non-const access to dictionary
        dictionary& dict();

        //- Write
        void write(Ostream&) const;

        //- Return info proxy.
        //  Used to print token information to a stream
        InfoProxy<dictionaryEntry> info() const
        {
            return *this;
        }


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const dictionaryEntry&);
};


template<>
Ostream& operator<<(Ostream&, const InfoProxy<dictionaryEntry>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
