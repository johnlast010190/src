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
    (c) 2012 OpenFOAM Foundation

Class
    Foam::simpleRegIOobject

Description
    Abstract base class for registered object with I/O. Used in debug symbol
    registration.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef simpleRegIOobject_H
#define simpleRegIOobject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                         Class simpleRegIOobject Declaration
\*---------------------------------------------------------------------------*/

class simpleRegIOobject
{
public:

    // Constructors

        //- Construct from objectregistry inserter and name
        simpleRegIOobject
        (
            void (*fn)(const char* name, simpleRegIOobject*),
            const char* name
        )
        {
            (*fn)(name, this);
        }


    //- Destructor
    virtual ~simpleRegIOobject()
    {};


    // Member Functions

        //- Read
        virtual void readData(Istream&) = 0;

        //- Write
        virtual void writeData(Ostream&) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
