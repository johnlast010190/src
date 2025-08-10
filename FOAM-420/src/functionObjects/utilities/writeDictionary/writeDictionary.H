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
    (c) 2016 OpenCFD Ltd.
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::functionObjects::writeDictionary

Group
    grpUtilitiesFunctionObjects

Description
    Writes dictionaries on start-up and on change.

SourceFiles
    writeDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeDictionary_H
#define functionObjects_writeDictionary_H

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"
#include "primitives/strings/lists/wordList.H"
#include "primitives/hashes/SHA1/SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeDictionary Declaration
\*---------------------------------------------------------------------------*/

class writeDictionary
:
    public regionFunctionObject
{
    // Private data

        //- Names of dictionaries to monitor
        wordList dictNames_;

        //- List of changed dictionaries (only those registered to database)
        List<SHA1Digest> digests_;


    // Private Member Functions

        //- Helper function to write the dictionary if found at location
        bool tryDirectory
        (
            const label dicti,
            const word& location,
            bool& firstDict
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        writeDictionary(const writeDictionary&);

        //- Disallow default bitwise assignment
        void operator=(const writeDictionary&);


public:

    //- Runtime type information
    TypeName("writeDictionary");


    // Constructors

        //- Construct from Time and dictionary
        writeDictionary
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeDictionary();


    // Member Functions

        //- Read the writeDictionary data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the selected dictionaries
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
