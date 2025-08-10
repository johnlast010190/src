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
    Foam::radiation::multiBandSolidAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description

    multiBandSolidAbsorptionEmission radiation absorption/emission for solids.


SourceFiles
    multiBandSolidAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef multiBandSolidAbsorptionEmission_H
#define multiBandSolidAbsorptionEmission_H

#include "submodels/absorptionEmissionModel/absorptionEmissionModel/absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class multiBandSolidAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class multiBandSolidAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of bands
        static const int maxBands_ = 5;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //-  Absorption coefficients
        scalarList absCoeffs_;

        //-  Emissivity coefficients
        scalarList emiCoeffs_;

        //- Bands
        label nBands_;


public:

    //- Runtime type information
    TypeName("multiBandSolidAbsorptionEmission");


    // Constructors

        //- Construct from components
        multiBandSolidAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~multiBandSolidAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient
                tmp<volScalarField> aCont(const label bandI) const;


            // Emission coefficient

                //- Emission coefficient
                tmp<volScalarField> eCont(const label bandI) const;


            // Emission contribution

                //- Emission contribution
                tmp<volScalarField> ECont(const label bandI) const;


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
