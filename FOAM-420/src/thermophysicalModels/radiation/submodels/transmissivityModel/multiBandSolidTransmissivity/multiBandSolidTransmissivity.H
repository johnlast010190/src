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
    (c) 2015 OpenCFD Ltd.

Class
    Foam::radiation::multiBandSolidTransmissivity

Group
    grpRadiationTransmissivitySubModels

Description
    multiBandSolidTransmissivity radiation transmissivity for solids.

SourceFiles
    multiBandSolidTransmissivity.C

\*---------------------------------------------------------------------------*/

#ifndef multiBandSolidTransmissivity_H
#define multiBandSolidTransmissivity_H

#include "submodels/transmissivityModel/transmissivityModel/transmissivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class multiBandSolidTransmissivity Declaration
\*---------------------------------------------------------------------------*/

class multiBandSolidTransmissivity
:
    public transmissivityModel
{
public:

    // Public data

        //- Maximum number of bands
        static const label maxBands_ = 5;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //-  Absorption coefficients
        scalarList tauCoeffs_;

        //- Bands
        label nBands_;


public:

    //- Runtime type information
    TypeName("multiBandSolidTransmissivity");


    // Constructors

        //- Construct from components
        multiBandSolidTransmissivity
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~multiBandSolidTransmissivity();


    // Member Functions

        // Access

            //- Transmissivity coefficient
            tmp<volScalarField> tauEff(const label bandI) const;


            //- Is Grey
            inline bool isGrey() const
            {
                return false;
            }


            //- Number of bands
            inline label nBands() const
            {
                return nBands_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
