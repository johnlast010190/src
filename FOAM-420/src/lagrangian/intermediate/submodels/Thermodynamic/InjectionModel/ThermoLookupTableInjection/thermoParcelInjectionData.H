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
    Foam::thermoParcelInjectionData

Description
    Container class to provide injection data for thermodynamic parcels

SourceFiles
    thermoParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef thermoParcelInjectionData_H
#define thermoParcelInjectionData_H

#include "submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class thermoParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const thermoParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    thermoParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
                 Class thermoParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class thermoParcelInjectionData
:
    public kinematicParcelInjectionData
{
protected:

    // Parcel properties

        //- Temperature [K]
        scalar T_;

        //- Specific heat capacity [J/kg/K]
        scalar Cp_;


public:

    //- Runtime type information
    TypeName("thermoParcelInjectionData");

    // Constructors

        //- Null constructor
        thermoParcelInjectionData();

        //- Construct from dictionary
        thermoParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        thermoParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~thermoParcelInjectionData();


    // Access

        //- Return const access to the temperature
        inline scalar T() const;

        //- Return const access to the specific heat capacity
        inline scalar Cp() const;


    // Edit

        //- Return access to the temperature
        inline scalar& T();

        //- Return access to the specific heat capacity
        inline scalar& Cp();


    // I-O

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const thermoParcelInjectionData& data
        );

        //- Istream operator
        friend Istream& operator>>
        (
            Istream& is,
            thermoParcelInjectionData& data
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
