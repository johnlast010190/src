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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Class
    Foam::functionObjects::runTimeControls::runTimeCondition

Description
    Base class for run time conditions

SourceFiles
    runTimeCondition.C
    runTimeConditionNew.C
    runTimeCondition.H

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeControls_runTimeCondition_H
#define functionObjects_runTimeControls_runTimeCondition_H

#include "db/functionObjects/stateFunctionObject/stateFunctionObject.H"
#include "db/dictionary/dictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                      Class runTimeCondition Declaration
\*---------------------------------------------------------------------------*/

class runTimeCondition
{

protected:

    // Protected data

        //- Condition name
        word name_;

        //- Reference to the object registry
        const objectRegistry& obr_;

        //- State
        stateFunctionObject& state_;

        //- On/off switch
        bool active_;

        //- Reference to the condition dictionary
        dictionary& conditionDict_;

        //- Group index - if applied, all conditions in a group must be
        //  satisfield before condition is met
        label groupID_;

public:
    // Public data

        //- Switch to send output to Info
        Switch& log;


protected:

    // Protected Member Functions

        //- Set the condition dictionary (create if necessary)
        dictionary& setConditionDict();

        //- Return const access to the conditions dictionary
        const dictionary& conditionDict() const;

        //- Return non-const access to the conditions dictionary
        dictionary& conditionDict();


public:

    //- Runtime type information
    TypeName("runTimeCondition");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        runTimeCondition,
        dictionary,
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            stateFunctionObject& state
        ),
        (name, obr, dict, state)
    );


    //- Constructor
    runTimeCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~runTimeCondition();

    //- Selector
    static autoPtr<runTimeCondition> New
    (
        const word& conditionName,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );


    // Public Member Functions

        //- Return the condition name
        virtual const word& name() const;

        //- Return the active flag
        virtual bool active() const;

        //- Return the group index
        virtual label groupID() const;

        //- Apply the condition
        virtual bool apply() = 0;

        //- Helper functions to get/store values in the state specific to a condition

            template <class Type>
            void setConditionResult
            (
                const word& entryName,
                const Type& value
            );

            template <class Type>
            Type getConditionResult
            (
                const word& entryName,
                const Type& defaultValue
            ) const;

            template <class Type>
            Type getConditionResult
            (
                const word& entryName
            ) const;


        //- Write
        virtual void write() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "runTimeControl/runTimeCondition/runTimeCondition/runTimeConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
