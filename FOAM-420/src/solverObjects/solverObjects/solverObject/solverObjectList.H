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

Class
    Foam::solverObjectList

Description
    Convenience class to gather together all the solverObjects in a region.

SourceFiles
    solverObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef solverObjectList_H
#define solverObjectList_H

#include "solverObject.H"
#include "solverOption/solverOption.H"
#include "containers/Lists/UPtrList/UPtrList.H"
#include "cfdTools/general/fvOptions/fvOption.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solverObjectList Declaration
\*---------------------------------------------------------------------------*/

class solverObjectList
:
    public UPtrList<solverObject>
{
protected:

    //- The solver options from which the solver objects originated
    UPtrList<fv::solverOption> solverOptions_;

    //- List of correctors for each solver object
    List<HashSet<word>> correctorNames_;

    fv::optionList& fvOptions_;

public:

    //- Runtime type information
    TypeName("solverObjectList");

    // Constructors

        //- Construct by pulling out fv::option pointers from fvOptionList
        solverObjectList(fv::optionList& fvOptions, bool legacySolver);

    // Destructor

        virtual ~solverObjectList()
        {}

    // Member Functions

        virtual void getSolveGraph
        (
            DynamicList<word>& solveNames,
            HashTable<wordList>& derivedFields,
            SolveTable<List<Tuple2<solveID, bool>>>& solveDependencies,
            DynamicList<label>& solverOrder,
            HashTable<wordList>& correctorMembers
        );

        virtual scalar getMaxTimeStep();

        virtual bool isFinalCorrector
        (
            const label corrector,
            const word& correctorName
        );

        void beginIteration
        (
            const label corrector,
            const word& correctorName,
            const bool finalIter
        );

        virtual void endIteration
        (
            const label corrector,
            const word& correctorName,
            const bool finalIter
        );

        virtual void updateTimings();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
