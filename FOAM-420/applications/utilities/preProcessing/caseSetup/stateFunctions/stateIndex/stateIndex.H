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
    (c) 2010-2020 Esi Ltd
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::stateIndex

Description
    Determines the type of solver based on dictionary settings and
    ouputs appropriate dictionaries.
    Also contains field prototypes and associated processing

SourceFiles
    stateIndex.C

\*---------------------------------------------------------------------------*/

#ifndef stateIndex_H
#define stateIndex_H

#include "containers/Lists/DynamicList/DynamicList.H"
#include "stateFunction/stateFunction.H"
#include "stateFunctions/global/globalState.H"
#include "solutionRegions/solutionScheduler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class stateIndex Declaration
\*---------------------------------------------------------------------------*/

class stateIndex
{

private:

    // Private data

        // Global references & dictionary input data

            //- const reference to time
            const Time& runTime_;

            //- region list
            HashTable<wordList> regionGroups_;

            solutionScheduler solutionSchedule_;

            //- global input handling
            stateFunctions::globalState global_;

        //- regional input assembly

            //- regional state data
            PtrList<stateFunction> regions_;


    // Private Member Functions


        //- Extract region information and output regionProperties file
        void regionProperties();

        //- combine explicit, group and regExp input
        Xfer<dictionary> assembleRegionInput
        (
            const dictionary& input,
            const word& regionName,
            const word& rgID
        ) const;


        //- compile input ditionaries
        void assembleDictionaries
        (
            const dictionary& input,
            dictionary& defaults
        );

        //- load libraries from in-state controlDict
        void openLibs();


        //- Disallow default bitwise copy construct and assignment
        stateIndex(const stateIndex&);
        void operator=(const stateIndex&);


public:

    // Constructors

        //- Construct from dicitonary
        stateIndex
        (
            const Time& runTime,
            const dictionary& input,
            dictionary& defaults,
            const bool distributed,
            const bool collated
        );


    // Member Functions


        // Edit

            //- write dictionaries to file
            void writeDictionaries();

            //- create mesh and store old fields
            void createRegionObjects();

            //- modify mesh
            void modifyMesh(scalar writePause);

            //- initialise fields
            void createFields();

            //- switch reset boundaries
            void resetBoundaries(const bool& reset);

        // Access

            const PtrList<stateFunction>& regions() const
            {
                return regions_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
