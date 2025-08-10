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
    (c) 2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::processorRunTimes

Description

SourceFiles
    processorRunTimes.C

\*---------------------------------------------------------------------------*/

#ifndef processorRunTimes_H
#define processorRunTimes_H

#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class processorRunTimes Declaration
\*---------------------------------------------------------------------------*/

class processorRunTimes
{
    // Private Data

        //- The complete run time
        Time completeRunTime_;

        //- Optional alternative decomposeParDict location
        const fileName decomposeDictFile_;

        //- Number of processors
        const label nProcs_;

        //- Processor run times
        PtrList<Time> procRunTimes_;


public:

    // Constructors

        //- Construct from directory and arguments (used in reconstructPar)
        processorRunTimes(const word& name, const argList& args);

        //- Construct from directory, arguments and optional alternative
        //  decomposition dictionary location (used in decomposePar)
        processorRunTimes
        (
            const word& name,
            const argList& args,
            const fileName& decomposeDictFile
        );


    //- Destructor
    ~processorRunTimes();


    // Member Functions

        //- Access the complete run time
        const Time& completeTime() const
        {
            return completeRunTime_;
        }

        //- Access the decomposeParDict location
        const fileName& decomposeDictFile() const
        {
            return decomposeDictFile_;
        }

        //- Return the number of processors
        label nProcs() const
        {
            return nProcs_;
        }

        //- Access the processor run times
        const PtrList<Time>& procTimes() const
        {
            return procRunTimes_;
        }

        //- Set the time
        void setTime(const instant& inst, const label newIndex);

        //- Select the time
        instantList selectComplete(const argList& args);

        //- Select the time
        instantList selectProc(const argList& args);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //