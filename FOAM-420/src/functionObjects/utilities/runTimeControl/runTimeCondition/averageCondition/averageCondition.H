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

Class
    Foam::functionObjects::runTimeControls::averageCondition

Description
    Average run time condition - satisfied when average does not change by
    more than a given value.

SourceFiles
    averageCondition.H
    averageCondition.C

\*---------------------------------------------------------------------------*/

#ifndef averageCondition_H
#define averageCondition_H

#include "runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "primitives/bools/Switch/Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                      Class averageCondition Declaration
\*---------------------------------------------------------------------------*/

class averageCondition
:
    public runTimeCondition,
    public writeFile
{
protected:

    // Protected data

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of fields on which to operate
        wordList fieldNames_;

        //- Satisfied when difference in mean values is less than this value
        const scalar tolerance_;

        //- Averaging window
        const scalar window_;

        //- StartTime discard
        const scalar startTime_;

        //- Average time per field
        List<scalar> totalTime_;

        //- Reset the averaging process on restart flag
        Switch resetOnRestart_;


    // Protected Member Functions

        //- Templated function to calculate the average
        template<class Type>
        void calc
        (
            const word& fieldName,
            const scalar alpha,
            const scalar beta,
            bool& satisfied,
            bool& processed
        );


public:

    //- Runtime type information
    TypeName("average");

    //- Constructor
    averageCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~averageCondition();


    // Public Member Functions

        //- Apply the condition
        virtual bool apply();

        //- Write
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "runTimeControl/runTimeCondition/averageCondition/averageConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
