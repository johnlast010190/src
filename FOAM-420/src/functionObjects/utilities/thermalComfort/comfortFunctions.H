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
    (c) 2010-2011, Esi Ltd

Class
    Foam::comfortFunctions

Description
    A set of functions for evaluating comfort

SourceFiles
    comfortFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef comfortFunctions_H
#define comfortFunctions_H

#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class comfortFunctions Declaration
\*---------------------------------------------------------------------------*/

class comfortFunctions
{
    // Private data

        //- global resistivity clothes
        const scalar clo_;

        //- global metabolic rate
        const scalar met_;

        //- global work
        const scalar work_;

        //Steffen-Boltzman contant
        const scalar sigma_;

    public:

        ClassName("comfortFunctions");

        //Constructor used for ISO7730 and ASHRAE55 calculation
        comfortFunctions
        (
            const scalar clo,
            const scalar met,
            const scalar work,
            const scalar Tu
        );

        //- Wet bulb temperature from inverse solution method of Roland Stull
        scalar wetBulb
        (
            const scalar& temperature,
            const scalar& relativeHumidity
        );

        //- Return Saturation Vapour Pressure (Pa)
        scalar svp
        (
            const scalar Ta // K
        );

        //- Return Saturation Vapour Pressure (Torr)
        scalar svpTorr
        (
            const scalar Ta // K
        );

        //- Calculation of Operative Temperature
        scalar TOp
        (
            const scalar& va, // ms-1
            const scalar& temperature, // K
            const scalar& radiantTemperature // K
        );

        //- Calculation of Predicted Percentage Dissatisfied (PPD)
        scalar PPD
        (
            const scalar& PMV
        );

        //- Calculation of Predicted Mean Vote (PMV)
        scalar PMV
        (
            const scalar& va, // ms-1
            const scalar& temperature, // K
            const scalar& radiantTemperature, // K
            const scalar& relativeHumidity // %
        );

        //- Calculation of Draft Rate
        scalar DR
        (
            const scalar& va,
            const scalar& temperature,
            const scalar& Tu_
        );

        scalar TU
        (
            const scalar& va,
            const scalar& k,
            const scalar upperLimitTu,
            const scalar lowerLimitTu
        );

        //- Calculation of Universal Thermal Comfort Index (UTCI)
        scalar UTCI
        (
            const scalar& va, // ms-1
            const scalar& temperature,// K
            const scalar& radiantTemperature,// K
            const scalar& relativeHumidity // %
        );

        //- Calculation of Standard Effective Temperatue (SET)
        scalar SET
        (
            const scalar& velocity, // ms-1
            const scalar& temperature, // K
            const scalar& pressure, // Nm-2
            const scalar& radiantTemperature,// K
            const scalar& relativeHumidity // %
        );

        //- Calculation of corrected Standard Effective Temperatue (SET)
        scalar correctSET
        (
            const scalar& pressure, // Nm-2
            const scalar& relativeHumidity, //%
            scalar& velocity, // ms-1
            scalar& temperature, // K
            scalar& radiantTemperature, // K
            const bool occupantControl
        );

        //- Calculation of PHS Model - ISO7933 fields
        scalarField PHSiso7933
        (
            const scalar& va,
            const scalar& temperature,
            const scalar& radiantTemperature,
            const scalar& relativeHumidity,
            const scalar& weight,
            const scalar& height,
            const Switch& drink,
            const scalar& duration,
            const scalar& posture,
            const Switch& accl,
            const scalar& thetaW,
            const scalar& walkSpeed
        );

        //- Tnwb calculation - WBGT Model - ISO7243 fields
        scalar fTnwb
        (
            scalar& Tit,
            const scalar& pressure,
            const scalar& temperature,
            const scalar& relativeHumidity,
            const scalar& va,
            const scalar& rho,
            const scalar& nu,
            const scalar& Cp,
            const scalar& kappa,
            const scalar& solarRad,
            const scalar& zenithAngle
        );

        //- Tglobe calculation - WBGT Model - ISO7243 fields
        //- returns the residual from energy balance equation applied at
        //- black globe temperature probe
        scalar fTg
        (
            scalar& Tit,
            const scalar& temperature,
            const scalar& relativeHumidity,
            const scalar& va,
            const scalar& rho,
            const scalar& nu,
            const scalar& Cp,
            const scalar& kappa,
            const scalar& solarRad,
            const scalar& zenithAngle
        );

        //- Calculation of WBGT Model - ISO7243 fields
        scalar WBGTiso7243
        (
            const scalar& pressure,
            const scalar& va,
            const scalar& temperature,
            const scalar& relativeHumidity,
            const scalar& wetBulbTemperature,
            const scalar& rho,
            const scalar& nu,
            const scalar& Cp,
            const scalar& kappa,
            const scalar& solarRad,
            const Switch& flagSolarRad,
            const scalar& zenithAngle
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
