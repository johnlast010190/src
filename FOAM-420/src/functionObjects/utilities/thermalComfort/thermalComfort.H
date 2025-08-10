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
    (c) 2013 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::functionObjects::thermalComfort

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates thermal comfort fields

    Wet bulb temperature from inverse solution method of Roland Stull
    valid for relative humidity range [5,99]% and temperature trange [-20,50]C

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    thermalComfort.C
    IOthermalComfort.H

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_thermalComfort_H
#define functionObjects_thermalComfort_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class thermalComfort Declaration
\*---------------------------------------------------------------------------*/

class thermalComfort
:
    public fvMeshFunctionObject
{
    protected:
        //- Return Cp
        virtual tmp<volScalarField> getCp() const;

        //- Return nu - kinematic viscosity
        virtual tmp<volScalarField> getNu() const;

        //- Return kappa - thermal conductivity
        virtual tmp<volScalarField> getKappa() const;

        //- Return Rho
        virtual tmp<volScalarField> getRho() const;

        //- Construction of fields
        void storeField
        (
            const word& fieldName,
            const dimensionSet& dims = dimless,
            const scalar initValue = Zero
        );

    private:
        // Private data

        //- Whether to calculate UTCI field
        Switch calculateUTCI_;

        //- Whether to correct local velocity to 10m height
        Switch correctUTCI_;

        //- Whether to calculate ISO 7730 fields PPD, PPE, DR
        Switch calculateISO7730_;

        //- Whether to calculate AHRAE-55-17 fields
        Switch calculateASHRAE55_;

        //- Whether to calculate Standard Effective Temperature (SET)
        Switch calculateASHRAE55PMV_;

        //- Whether to calculate wet bulb temperature
        Switch calculateWetBulbTemp_;

        //- Whether to calculate ISO 7933 fields - PHS Model
        Switch calculateISO7933_;

        //- Whether to calculate ISO 7243 fields - WBGT Model
        Switch calculateISO7243_;

        //- resistivity clothes
        scalar clo_;

        //- metabolic rate
        scalar met_;

        //- work
        scalar work_;

        //- average weight of human body
        scalar weight_;

        //- average height of human body
        scalar height_;

        //- Whether workes can drink water as they want or not
        Switch drink_;

        //- Duration of the work sequence in minutes
        scalar duration_;

        //- posture = 1 (sitting) / posture = 2 (standing) / posture = 3 (crouching)
        scalar posture_;

        //- Whether the workers are acclimatised or not (used to the environment)
        Switch accl_;

        //- angle between walking direction and wind direction [deg]
        scalar thetaW_;

        //- walking speed
        scalar walkSpeed_;

        //- Whether occupant controlled room (ASHRAE-55-17 standard)
        Switch occupantControl_;

        //- Relative humidity (if no water vapour field)
        scalar RH_;

        //- Turbulent intensity
        scalar Tu_;

        //- Whether to calculate Turbulent Intensity (Tu) locally
        Switch calculateTu_;

        //- Upper Limit of Turbulent Intensity (Tu) when calculated locally
        scalar upperLimitTu_;

        //- Lower Limit of Turbulent Intensity (Tu) when calculated locally
        scalar lowerLimitTu_;

        //- Name of velocity field, default is "U"
        word UName_;

        //- Name of pressure field, default is "p"
        word pName_;

        //- List of solar fields
        wordList solarFields_;

        //- Reference pressure
        autoPtr<dimensionedScalar> pRef_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalComfort(const thermalComfort&);

        //- Disallow default bitwise assignment
        void operator=(const thermalComfort&);


public:

    //- Runtime type information
    TypeName("thermalComfort");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        thermalComfort
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermalComfort();


    // Member Functions

        //- Return reference pressure
        dimensionedScalar pRef();

        //- Create thermal comfort fields
        void createFields();

        //- Read the thermalComfort data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate the thermalComfort and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
