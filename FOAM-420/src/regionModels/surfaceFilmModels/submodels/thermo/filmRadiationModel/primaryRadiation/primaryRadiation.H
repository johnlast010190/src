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
    (c) 2012-2017 OpenFOAM Foundation

Class
    Foam::regionModels::surfaceFilmModels::primaryRadiation

Description
    Radiation model whereby the radiative heat flux is mapped from the primary
    region

SourceFiles
    primaryRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef primaryRadiation_H
#define primaryRadiation_H

#include "submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModel.H"
#include "fields/volFields/volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class primaryRadiation Declaration
\*---------------------------------------------------------------------------*/

class primaryRadiation
:
    public filmRadiationModel
{
private:

    // Private data

        //- Incident radiative flux mapped from  the primary region / [kg/s3]
        volScalarField qinPrimary_;


    // Private member functions

        //- Disallow default bitwise copy construct
        primaryRadiation(const primaryRadiation&);

        //- Disallow default bitwise assignment
        void operator=(const primaryRadiation&);


public:

    //- Runtime type information
    TypeName("primaryRadiation");


    // Constructors

        //- Construct from surface film model and dictionary
        primaryRadiation
        (
            surfaceFilmModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~primaryRadiation();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct();

            //- Return the radiation sensible enthalpy source
            //  Also updates qrNet
            virtual tmp<volScalarField> Shs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
