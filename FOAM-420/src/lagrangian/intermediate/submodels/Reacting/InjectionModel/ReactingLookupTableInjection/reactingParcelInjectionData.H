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
    Foam::reactingParcelInjectionData

Description
    Container class to provide injection data for reacting parcels

SourceFiles
    reactingParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef reactingParcelInjectionData_H
#define reactingParcelInjectionData_H

#include "submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionData.H"
#include "primitives/Scalar/lists/scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class reactingParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const reactingParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    reactingParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
               Class reactingParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class reactingParcelInjectionData
:
    public thermoParcelInjectionData
{
protected:

    // Parcel properties

        //- List of mass fractions
        scalarList Y_;


public:

    //- Runtime type information
    TypeName("reactingParcelInjectionData");

    // Constructors

        //- Null constructor
        reactingParcelInjectionData();

        //- Construct from dictionary
        reactingParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        reactingParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~reactingParcelInjectionData();


    // Access

        //- Return const access to the list of mass fractions
        inline const scalarList& Y() const;


    // Edit

        //- Return access to the mass fractions
        inline scalarList& Y();


    // I-O

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const reactingParcelInjectionData& data
        );

        //- Istream operator
        friend Istream& operator>>
        (
            Istream& is,
            reactingParcelInjectionData& data
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "submodels/Reacting/InjectionModel/ReactingLookupTableInjection/reactingParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
