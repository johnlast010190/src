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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::reactionRateFlameAreaModels::relaxation

Description
    Consumption rate per unit of flame area obtained from a relaxation equation

SourceFiles
    relaxation.C

\*---------------------------------------------------------------------------*/

#ifndef relaxation_H
#define relaxation_H

#include "FSD/reactionRateFlameAreaModels/reactionRateFlameArea/reactionRateFlameArea.H"
#include "FSD/reactionRateFlameAreaModels/consumptionSpeed/consumptionSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRateFlameAreaModels
{

/*---------------------------------------------------------------------------*\
                           Class relaxation Declaration
\*---------------------------------------------------------------------------*/

class relaxation
:
    public reactionRateFlameArea
{
    // Private Data

        //- Correlation
        consumptionSpeed correlation_;

        //- Proportionality constant for time scale in the relaxation Eq.
        scalar C_;

        //- Proportionality constant for sub-grid strain
        scalar alpha_;


     // Private Member Functions

        //- Disallow copy construct
        relaxation(const relaxation&);

        //- Disallow default bitwise assignment
        void operator=(const relaxation&);


public:

    //- Runtime type information
    TypeName("relaxation");


    // Constructors

        //- Construct from dictionary and psiReactionThermo
        relaxation
        (
            const word modelType,
            const dictionary& dictCoeffs,
            const objectRegistry& obr,
            const combustionModel& combModel
        );


    // Destructor

        virtual ~relaxation();


    // Member functions

        //- Correct omega
        virtual void correct(const volScalarField& sigma);


    // IO

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dictProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End reactionRateFlameAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
