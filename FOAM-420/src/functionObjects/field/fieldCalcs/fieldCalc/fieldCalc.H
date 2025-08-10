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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Class
    Foam::functionObjects::fieldCalc

Group
    grpFieldFunctionObjects

Description
    Base class for field value-based function objects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    fieldCalc.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldCalc_H
#define functionObjects_fieldCalc_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "primitives/bools/Switch/Switch.H"
#include "fields/Fields/Field/Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldCalc Declaration
\*---------------------------------------------------------------------------*/

class fieldCalc
:
    public fvMeshFunctionObject
{

protected:

    // Protected data

        //- Construction dictionary
        dictionary dict_;

        //- Name of region (patch, zone, etc.)
        word regionName_;

    // Protected Member Functions

        //- Combine fields from all processor domains into single field
        template<class Type>
        void combineFields(Field<Type>& field);

        //- Combine fields from all processor domains into single field
        template<class Type>
        void combineFields(tmp<Field<Type>>&);


public:

    //- Run-time type information
    TypeName("fieldCalc");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fieldCalc,
            dictionary,
            (
                const word& name,
                const objectRegistry& obr,
                const dictionary& dict
            ),
            (name, obr, dict)
        );


    // Constructors

        //- Construct from Time and dictionary
        fieldCalc
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word& valueType
        );

        //- Construct from objectRegistry and dictionary
        fieldCalc
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const word& valueType
        );

        //- Return a reference to the selected fieldCalc
        static autoPtr<fieldCalc> New
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool output = true
        );


    //- Destructor
    virtual ~fieldCalc();


    // Member Functions

        //- Return the reference to the construction dictionary
        inline const dictionary& dict() const;

        //- Return the region name
        inline const word& regionName() const;

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldCalcs/fieldCalc/fieldCalcI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldCalcs/fieldCalc/fieldCalcTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
