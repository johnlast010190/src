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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::functionObjects::aniCourantNo

Group
    grpFieldFunctionObjects

Description
    This function object calculates and outputs the Courant number as a
    volScalarField.  The field is stored on the mesh database so that it can
    be retrieved and used for other applications.

Usage
    Example of function object specification to calculate the Courant number:
    \verbatim
    aniCourantNo1
    {
        type        aniCourantNo;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: aniCourantNo    | yes         |
        rho          | Name of density field   | no          | rho
        field        | Name of flux field      | no          | phi
        result       | Name of Courant number field | no     | <function name>
        log          | Log to standard output  | no          | yes
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    aniCourantNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_aniCourantNo_H
#define functionObjects_aniCourantNo_H

#include "fieldExpression/fieldExpression.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class aniCourantNo Declaration
\*---------------------------------------------------------------------------*/

class aniCourantNo
:
    public fieldExpression
{
    // Private data

        //- Name of density field (optional)
        word rhoName_;


        word UName_;

    // Private Member Functions

        //- Divide the Courant number by rho if required
        tmp<volScalarField::Internal> byRho
        (
            const tmp<volScalarField::Internal>& Co
        ) const;

        //- Calculate the Courant number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("aniCourantNo");


    // Constructors

        //- Construct from Time and dictionary
        aniCourantNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~aniCourantNo();


    // Member Functions

        //- Read the aniCourantNo data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
