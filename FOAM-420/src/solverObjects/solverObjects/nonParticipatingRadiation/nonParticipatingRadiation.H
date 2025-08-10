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
    (c) 2010-2019 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::nonParticipatingRadiation

Description
    A function object for solving non-participating radiation transport.
    Reads input from constant/radiationProperties

    Example of the function object specification:
    @verbatim
    radiation
    {
        type            nonParticipatingRadiation;
        region          region0;

        //T T;
    }
    @endverbatim

SourceFiles
    nonParticipatingRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef nonParticipatingRadiation_H
#define nonParticipatingRadiation_H

#include "solverObject/solverObject.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "cfdTools/general/include/fvCFD.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class nonParticipatingRadiation Declaration
\*---------------------------------------------------------------------------*/

class nonParticipatingRadiation
:
    public solverObject
{
protected:

    // Protected data

        //- temperature name
        word Tname_;

        //- radiation model
        autoPtr<radiation::radiationModel> radiation_;

    // Protected functions

        //- Disallow default bitwise copy construct
        nonParticipatingRadiation(const nonParticipatingRadiation&);

        //- Disallow default bitwise assignment
        void operator=(const nonParticipatingRadiation&);


public:

    //- Runtime type information
    TypeName("nonParticipatingRadiation");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        nonParticipatingRadiation
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );


    // Destructor

        virtual ~nonParticipatingRadiation()
        {}


    // Member Functions

        //- Read the parameters
        virtual void read(const dictionary&);

        //- Identify any solves handled by this solverObject, their levels of
        // iteration and dependent solves
        virtual void getSolveGraph
        (
            DynamicList<word>& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        )
        {
            solveNames.append("I");
            correctorMembers.insert
            (
                solverObject::outerCorrectorName, solveNames
            );
        }

        //- Update solver
        virtual void correct(const word& solveName);

        //- Sample and write
        virtual void write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
