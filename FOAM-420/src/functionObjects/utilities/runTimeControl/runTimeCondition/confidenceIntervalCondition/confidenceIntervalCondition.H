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
    (c) 2018 Esi Ltd.
    (c) 2015 OpenFOAM Foundation

Class
    Foam::functionObjects::runTimeControls::confidenceIntervalCondition

Description
    Confidence interval run time condition
    Satisfied when confidence interval drops below a specified range.

SourceFiles
    confidenceIntervalCondition.H
    confidenceIntervalCondition.C

\*---------------------------------------------------------------------------*/

#ifndef confidenceIntervalCondition_H
#define confidenceIntervalCondition_H

#include "runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "primitives/enums/NamedEnum.H"
#include "primitives/bools/Switch/Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                      Class confidenceIntervalCondition Declaration
\*---------------------------------------------------------------------------*/

class confidenceIntervalCondition
:
    public runTimeCondition,
    public writeFile
{
public:

        //- Interval type options
        enum intvTypes
        {
            itRel,
            itAbs
        };


protected:

    // Protected data

        //- Interval type names
        static const NamedEnum<intvTypes,2> intvTypeNames_;

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of entries on which to operate
        wordList entryNames_;

        //- Required confidence - 0.01 -- 0.999
        const scalar confidence_;

        //- Termination interval
        const scalar interval_;

        //- interval mode: fraction/absolute
        const intvTypes intervalType_;

        //- Minimum samples
        label minSamples_;

        //- Minimum sample time
        scalar minSampleTime_;

        //- table file name
        fileName tDistFile_;

        //- Student-T distribution coefficient
        scalar Cst_;


    // Protected Member Functions

        //- function to calculate scalar statistics
        void calc
        (
            const word& entryName,
            const scalar& dt,
            bool& satisfied,
            bool& processed
        );

        //- lookup Student-T distribution coefficient
        scalar lookupSTD(const scalar& conf) const;


public:

    //- Runtime type information
    TypeName("confidenceInterval");

    //- Constructor
    confidenceIntervalCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~confidenceIntervalCondition();


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
    #include "runTimeControl/runTimeCondition/confidenceIntervalCondition/confidenceIntervalConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
