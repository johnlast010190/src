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
    Foam::unwatchedIOdictionary

Description
    unwatchedIOdictionary is like IOdictionary but stores
    dependencies as files instead of fileMonitor watchIndices. Used
    to read controlDict since there fileMonitor not yet setup.

\*---------------------------------------------------------------------------*/

#ifndef unwatchedIOdictionary_H
#define unwatchedIOdictionary_H

#include "db/IOobjects/IOdictionary/baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class unwatchedIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class unwatchedIOdictionary
:
    public baseIOdictionary
{

    fileNameList files_;

public:

    // Constructors

        //- Construct given an IOobject
        unwatchedIOdictionary(const IOobject& io);

        //- Construct given an IOobject and dictionary
        unwatchedIOdictionary
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Construct given an IOobject and Istream
        unwatchedIOdictionary(const IOobject& io, Istream& is);

    //- Destructor
    virtual ~unwatchedIOdictionary();


    // Member functions

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }

        //- Add file watch on object (READ_IF_MODIFIED)
        virtual void addWatch();

        //- Add file watch for fileName on object if not yet watched. Return
        //  index of watch
        virtual label addWatch(const fileName&);

        const fileNameList& files() const
        {
            return files_;
        }

        fileNameList& files()
        {
            return files_;
        }

};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<unwatchedIOdictionary>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
