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
    (c) 2018 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::baseIOdictionary

Description
    baseIOdictionary is derived from dictionary and IOobject to give the
    dictionary automatic IO functionality via the objectRegistry.

    To facilitate IO, baseIOdictionary is provided with a constructor from
    IOobject and with readData/writeData functions.

SourceFiles
    baseIOdictionary.C
    baseIOdictionaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef baseIOdictionary_H
#define baseIOdictionary_H

#include "db/dictionary/dictionary.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class baseIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class baseIOdictionary
:
    public regIOobject,
    public dictionary
{
    // Private data

        static bool writeDictionaries;

public:

    TypeName("dictionary");


    // Constructors

        //- Copy construct
        baseIOdictionary(const baseIOdictionary&) = default;

        //- Construct given an IOobject
        baseIOdictionary(const IOobject&);

        //- Construct given an IOobject and dictionary
        baseIOdictionary(const IOobject&, const dictionary&);

        //- Construct given an IOobject and Istream
        baseIOdictionary(const IOobject&, Istream&);


    //- Destructor
    virtual ~baseIOdictionary();


    // Member functions

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const = 0;

        //- Name function is needed to disambiguate those inherited
        //  from regIOobject and dictionary
        const word& name() const;

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const = 0;


    // Member operators

        //- Assignment of other baseIOdictionary's entries to this
        //  baseIOdictionary
        void operator=(const baseIOdictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
