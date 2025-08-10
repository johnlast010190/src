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
    (c) 2015-17 OpenFOAM Foundation

Class
    Foam::localIOdictionary

Description
    localIOdictionary is derived from IOdictionary but excludes parallel
    master reading.

SourceFiles
    localIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef localIOdictionary_H
#define localIOdictionary_H

#include "db/IOobjects/IOdictionary/baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class localIOdictionary
:
    public baseIOdictionary
{

public:

    // Constructors

        //- Construct given an IOobject
        localIOdictionary(const IOobject& io);

        //- Construct given an IOobject and dictionary
        localIOdictionary(const IOobject& io, const dictionary& dict);

        //- Construct given an IOobject and Istream
        localIOdictionary(const IOobject& io, Istream& is);

        //- Construct given an IOobject, supply wanted typeName
        localIOdictionary(const IOobject& io, const word& wantedType);


    //- Destructor
    virtual ~localIOdictionary();


    // Member functions

        //- Is object global
        virtual bool global() const
        {
            return false;
        }

        //- Return complete path + object name if the file exists
        //  in the case otherwise null
        virtual fileName filePath() const
        {
            // Use default (local only) search strategy
            return localFilePath(type());
        }

};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<localIOdictionary>()
{
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
