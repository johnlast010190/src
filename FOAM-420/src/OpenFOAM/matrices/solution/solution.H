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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2019-2022 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

Class
    Foam::solution

Description
    Selector class for relaxation factors, solver type and solution.

SourceFiles
    solution.C

\*---------------------------------------------------------------------------*/

#ifndef solution_H
#define solution_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class solution Declaration
\*---------------------------------------------------------------------------*/

class solution
:
    public IOdictionary
{
    // Private data

        //- Dictionary of temporary fields to cache
        dictionary cache_;

        //- Switch for the caching mechanism
        bool caching_;

        //- Dictionary of relaxation factors for all the fields
        dictionary fieldRelaxDict_;

        //- Cache of Function1s in above dictionary
        mutable HashPtrTable<Function1<scalar>> fieldRelaxCache_;

        //- Dictionary of relaxation factors for all the equations
        dictionary eqnRelaxDict_;

        //- Cache of Function1s in above dictionary
        mutable HashPtrTable<Function1<scalar>> eqnRelaxCache_;

        //- Optional default relaxation factor for all the fields
        autoPtr<Function1<scalar>> fieldRelaxDefault_;

        //- Optional default relaxation factor for all the equations
        autoPtr<Function1<scalar>> eqnRelaxDefault_;

        //- Dictionary of solver parameters for all the fields
        dictionary solvers_;

        //- Dictionary of direct set solution parameters
        dictionary localSolution_;

    // Private Member Functions

        //- Refresh solution dependent member data
        //- and apply local modifications
        void update();

public:
        //- Read settings from the dictionary
        void read(const dictionary&);

private:
        //- Disallow default bitwise copy construct and assignment
        solution(const solution&);
        void operator=(const solution&);


public:

    //- Update from older solver controls syntax
    //  Usually verbose, since we want to know about the changes
    //  Returns the number of settings changed
    static label upgradeSolverDict(dictionary& dict, const bool verbose=true);

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct for given objectRegistry and dictionary
        solution
        (
            const objectRegistry& obr,
            const fileName& dictName
        );


    // Member Functions

        // Access

            //- Return true if the given field should be cached
            bool cache(const word& name) const;

            //- Helper for printing cache message
            template<class FieldType>
            static void cachePrintMessage
            (
                const char* message,
                const word& name,
                const FieldType& vf
            );

            //- Return true if the relaxation factor is given for the field
            bool relaxField(const word& name) const;

            //- Return true if the relaxation factor is given and not one
            bool storeField(const word& name) const;

            //- Return true if the relaxation factor is given for the equation
            bool relaxEquation(const word& name) const;

            //- Return the relaxation factor for the given field
            scalar fieldRelaxationFactor(const word& name) const;

            //- Copy input to localSolution_
            // Call will force a reread of *this and apply
            // localSolution via a reverse merge operation to *this
            // after every read
            void setLocalSolutionDict(const dictionary&);


            //- Return const access to programatically fixed solution settings
            const dictionary& localSolutionDict() const
            {
                return localSolution_;
            }

            void resetReadOpt(const readOption& option);

            //- Return the relaxation factor for the given eqation
            scalar equationRelaxationFactor(const word& name) const;

            //- Return the selected sub-dictionary of solvers if the "select"
            //  keyword is given, otherwise return the complete dictionary
            const dictionary& dict() const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solverDict(const word& name) const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solver(const word& name) const;

            //- Check if solver dict exists
            bool solverDictExists(const word& name) const;


        // Read

            //- Read the solution dictionary
            bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/solution/solutionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
