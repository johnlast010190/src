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
    Foam::functionObjects::runTimeControls::equationMaxIterCondition

Description
    Maximum number of equation iterations run time condition

SourceFiles
    equationMaxIterCondition.H
    equationMaxIterCondition.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeConditions_equationMaxIterCondition_H
#define functionObjects_runTimeConditions_equationMaxIterCondition_H

#include "runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                  Class equationMaxIterCondition Declaration
\*---------------------------------------------------------------------------*/

class equationMaxIterCondition
:
    public runTimeCondition
{
protected:

    // Protected data

        //- Field name
        const wordList fieldNames_;

        //- Threshold for maximum number of iterations
        const label threshold_;

        //- Start checking from iteration - always skips first iteration
        label startIter_;


public:

    //- Runtime type information
    TypeName("equationMaxIter");

    //- Constructor
    equationMaxIterCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~equationMaxIterCondition();


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

#endif

// ************************************************************************* //
