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
    (c) 2019 Esi Ltd.

Class
    Foam::openPatch

Description
    Active patch if meanField on the patch is greater than a value.
    Used in applications with opening patches like doors/windows etc.

SourceFiles
    openPatch.C

\*---------------------------------------------------------------------------*/

#ifndef openPatch_H
#define openPatch_H

#include "patchSwitchFvMesh/patchSwitches/patchSwitch/patchSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class openPatch Declaration
\*---------------------------------------------------------------------------*/

class openPatch
:
    public patchSwitch
{
    // Private Member Functions

        //- Trigger alue that the GIB will get disabled if the mean value
        // of the field at the patch is fieldName_ is greater
        const scalar fieldLim_;

        //- Name of the field
        const word fieldName_;

        //- Side of the GIB (master/slave) that the condition is going to be
        //  evaluated
        const word side_;

        //- Evaluated patch ID
        label patchID_;

        //- The patch will never get enables.
        //  It is allowed to switch only from disable to enable
        virtual bool enableCondition()
        {
            return false;
        }

        virtual bool disableCondition();

        //- Disallow default bitwise copy construct
        openPatch(const openPatch&);

        //- Disallow default bitwise assignment
        void operator=(const openPatch&);


public:

    //- Runtime type information
    TypeName("openPatch");


    // Constructors

        //- Construct from objectRegistry, and read/write options
        openPatch
        (
            const fvMesh&,
            const dictionary&
        );

    // Selectors


    //- Destructor
    virtual ~openPatch(){};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
