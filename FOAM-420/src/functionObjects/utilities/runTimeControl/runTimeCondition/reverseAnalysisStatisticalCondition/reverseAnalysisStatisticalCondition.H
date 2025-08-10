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
    Foam::functionObjects::runTimeControls::reverseAnalysisStatisticalCondition

Description
    Confidence interval run time condition SOMETHING ELSE HERE.....
    Satisfied when confidence interval drops below a specified range.

SourceFiles
    reverseAnalysisStatisticalCondition.H
    reverseAnalysisStatisticalCondition.C

\*---------------------------------------------------------------------------*/

#ifndef reverseAnalysisStatisticalCondition_H
#define reverseAnalysisStatisticalCondition_H

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
                      Class reverseAnalysisStatisticalCondition Declaration
\*---------------------------------------------------------------------------*/

class reverseAnalysisStatisticalCondition
:
    public runTimeCondition,
    public writeFile
{
public:


protected:

    // Protected data

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of entries on which to operate
        wordList entryNames_;

        //- Minimum acceptable data points used in average
        const label minBox_;

        //- Standard deviation tolerance for asymptotic behaviour
        const scalar stdTol_;

        //- Tolerance on solution quality
        const scalar solQualTol_;

	//- Tolerance on gradient of best fit
        const scalar fitTol_;

        //- Settings used for discarding points
        const scalar discardRate_;

        //- Analyse vectors component wise
        const bool vectorComp_;

        //- Analyse vectors using their magnitude
        const bool vectorMag_;

        //- Analyse vectors using the sum of the components
        const bool vectorSum_;

        //- Starting flag for dictionary setting
        bool starting_;

	//- Tolerance on solution quality
	const bool resetOnRestart_;

    // Protected Member Functions

        //- function to calculate scalar statistics
	void calc
        (
            const word& entryName,
	    const scalar& currentValue,
            bool& satisfied,
            bool& processed
        );


public:

    //- Runtime type information
    TypeName("reverseAnalysisStatistical");

    //- Constructor
    reverseAnalysisStatisticalCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~reverseAnalysisStatisticalCondition();


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
    #include "runTimeControl/runTimeCondition/reverseAnalysisStatisticalCondition/reverseAnalysisStatisticalConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
