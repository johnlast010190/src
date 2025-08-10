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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Class
    Foam::functionObjects::magSqr

Group
    grpFieldFunctionObjects

Description
    Calculates the magnitude of the sqr of a field.

    The operation can be applied to any volume or surface field generating a
    volume or surface scalar field.
    With the %subRegion option, also supports fields on a surfMesh.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    magSqr.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_magSqr_H
#define functionObjects_magSqr_H

#include "fieldExpression/fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class magSqr Declaration
\*---------------------------------------------------------------------------*/

class magSqr
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the magnitude of the sqr of the field
        //  and register the result
        template<class Type>
        bool calcMagSqr();

        //- Calculate the magnitude of the sqr of the field
        //  and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("magSqr");


    // Constructors

        //- Construct from Time and dictionary
        magSqr
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~magSqr();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "magSqr/magSqrTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
