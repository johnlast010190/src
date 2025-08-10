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

Class
    Foam::Callback

Description
    Abstract class to handle automatic call-back registration with the
    CallbackRegistry.  Derive from this class and extend by adding the
    appropriate callback functions.

SourceFiles
    Callback.C

\*---------------------------------------------------------------------------*/

#ifndef Callback_H
#define Callback_H

#include "db/CallbackRegistry/CallbackRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Callback Declaration
\*---------------------------------------------------------------------------*/

template<class CallbackType>
class Callback
:
    public CallbackRegistry<CallbackType>::link
{
    // Private data

        //- CallbackRegistry reference
        CallbackRegistry<CallbackType>& cbr_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Callback<CallbackType>&);


public:

    // Constructors

        //- Construct and register with CallbackRegistry
        Callback(CallbackRegistry<CallbackType>&);

        //- Construct as copy
        Callback(const Callback<CallbackType>&);


    //- Destructor
    virtual ~Callback();


    // Member functions

        // Registration

            //- Register Callback with CallbackRegistry
            bool checkIn();

            //- Check-out Callback from CallbackRegistry
            bool checkOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "db/Callback/Callback.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
