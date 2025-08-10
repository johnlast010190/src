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
    (c) 2023 Esi Ltd.

Class
    Foam::fv::materialPropertiesSolver

Description
    Updates the material properties.

SourceFiles
    materialPropertiesSolver.C

\*---------------------------------------------------------------------------*/

#ifndef materialPropertiesSolver_H
#define materialPropertiesSolver_H

#include "solverObject/solverObject.H"
#include "solverOption/SolverOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class basicThermo;

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class materialPropertiesSolver Declaration
\*---------------------------------------------------------------------------*/

class materialPropertiesSolver
:
    public solverObject
{
    // Private data

        //- Thermo model
        basicThermo* thermoPtr_;


public:

    //- Runtime type information
    TypeName("materialPropertiesSolver");


    // Constructors

        //- Construct from components
        materialPropertiesSolver
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        materialPropertiesSolver(materialPropertiesSolver&) = delete;

    void operator=(materialPropertiesSolver&) = delete;


    // Solver Functions

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& derivedFields,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word&, const word&) override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
