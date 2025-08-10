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
    (c) 2015-2016 OpenCFD Ltd.

Class
    Foam::functionObjects::nullFunctionObject


Description
    Null function object for dealing with unavailable libraries
    Allows use of "exitIfNotFound false" option in functionObject
    run-time selection.

Usage
    \verbatim
    null
    {
        type            null;
        exitIfNotFound  false;
    }
    \endverbatim


    Adding "exitIfNotFound false;" to any function object entry will
    allow the selection to default to null function object.
    This can be useful during testing when some libraries might not be
    available or compiled.


SourceFiles
    nullfunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_nullFunctionObject_H
#define functionObjects_nullFunctionObject_H

#include "db/functionObjects/functionObject/functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class externalCoupled Declaration
\*---------------------------------------------------------------------------*/

class nullFunctionObject
:
    public functionObject
{

public:

    //- Runtime type information
    TypeName("null");

    // Constructors

        //- Construct given time and dictionary
        nullFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~nullFunctionObject();


    // Member Functions

      // Function object control

        //- Called at each ++ or += of the time-loop
        virtual bool execute();

        //- Write, currently a no-op
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
