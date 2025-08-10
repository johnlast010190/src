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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::solidThermo

Description
    Fundamental solid thermodynamic properties

SourceFiles
    solidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermo_H
#define solidThermo_H

#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fvMesh/fvMesh.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "basicThermo/basicThermo.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solidThermo Declaration
\*---------------------------------------------------------------------------*/

class solidThermo
:
    public basicThermo
{

protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'rhoThermo' to avoid (potential) conflict with solver density
        volScalarField rho_;


public:

    //- Runtime type information
    TypeName("solidThermo");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        objectRegistry,
        (const objectRegistry& obr, const word& phaseName),
        (obr, phaseName)
    );

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        dictionary,
        (const objectRegistry& obr, const dictionary& dict, const word& phaseName),
        (obr, dict, phaseName)
    );


    // Constructors

        //- Construct from obr and phase name
        solidThermo
        (
            const objectRegistry&,
            const word& phaseName
        );

        //- Construct from obr, dictionary and phase name
        solidThermo
        (
            const objectRegistry&,
            const dictionary& dict,
            const word& phaseName
        );

        //- Return a pointer to a new solidThermo created from
        //  the solidThermophysicalProperties dictionary
        static autoPtr<solidThermo> New
        (
            const objectRegistry&,
            const word& phaseName=word::null
        );

        //- Return a pointer to a new solidThermo created from
        //  local dictionary
        static autoPtr<solidThermo> New
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~solidThermo();


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Density [kg/m^3] (not applicable for solids)
            virtual tmp<volScalarField> buoyantRho() const
            {
                NotImplemented;
            }

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const = 0;

            //- Return true if thermal conductivity is isotropic
            virtual bool isotropic() const = 0;


        // Per patch calculation

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<vectorField> Kappa
            (
                const label patchi
            ) const = 0;


    // I-O

        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
