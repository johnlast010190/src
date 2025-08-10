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
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::regionModels::surfaceFilmModels::mappedConvectiveHeatTransfer

Description
    Convective heat transfer model based on a re-working of a Nusselt number
    correlation

SourceFiles
    mappedConvectiveHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef mappedConvectiveHeatTransfer_H
#define mappedConvectiveHeatTransfer_H

#include "submodels/thermo/heatTransferModel/heatTransferModel/heatTransferModel.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                Class mappedConvectiveHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class mappedConvectiveHeatTransfer
:
    public heatTransferModel
{
private:

    // Private data

        //- Heat transfer coefficient - primary region [W/m2/K]
        volScalarField htcConvPrimary_;

        //- Heat transfer coefficient - film region [W/m2/K]
        //  Assumes that the primary regtion to film region boundaries are
        //  described as mappedPushed types
        volScalarField htcConvFilm_;


    // Private member functions

        //- Disallow default bitwise copy construct
        mappedConvectiveHeatTransfer(const mappedConvectiveHeatTransfer&);

        //- Disallow default bitwise assignment
        void operator=(const mappedConvectiveHeatTransfer&);


public:

    //- Runtime type information
    TypeName("mappedConvectiveHeatTransfer");


    // Constructors

        //- Construct from surface film model and dictionary
        mappedConvectiveHeatTransfer
        (
            surfaceFilmModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~mappedConvectiveHeatTransfer();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct();

            //- Return the heat transfer coefficient [W/m2/K]
            virtual tmp<volScalarField> h() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
