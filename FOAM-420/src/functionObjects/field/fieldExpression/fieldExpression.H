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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::functionObjects::fieldExpression

Group
    grpFieldFunctionObjects

Description
    Base class for field expression function objects

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    fieldExpression.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldExpression_H
#define functionObjects_fieldExpression_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldExpression Declaration
\*---------------------------------------------------------------------------*/

class fieldExpression
:
    public fvMeshFunctionObject
{
protected:

    // Protected member data

        //- Name of field to process
        word fieldName_;

        //- Name of result field
        word resultName_;

        //- Optional constant value entry
        string valueStr_;

    // Protected member functions

        virtual bool calc() = 0;

        void setResultName(const word& typeName, const word& defaultArg);

        template<class Type>
        bool foundObject(const word& name, const bool verbose = true) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fieldExpression(const fieldExpression&);

        //- Disallow default bitwise assignment
        void operator=(const fieldExpression&);


public:

    //- Runtime type information
    TypeName("fieldExpression");


    // Constructors

        //- Construct from Time and dictionary
        fieldExpression
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word& fieldName = word::null,
            const word& resultName = word::null
        );


    //- Destructor
    virtual ~fieldExpression();


    // Member Functions

        //- Read the fieldExpression data
        virtual bool read(const dictionary&);

        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();

        //- Clear the result field from the objectRegistry
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldExpression/fieldExpressionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
