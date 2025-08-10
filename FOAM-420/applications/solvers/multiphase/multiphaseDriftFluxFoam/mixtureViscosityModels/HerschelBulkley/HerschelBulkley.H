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
    (c) 2021 Esi Ltd.

Class
    Foam::mixtureViscosityModels::HerschelBulkley

Description
     Viscosity correction model for HerschelBulkley.

SourceFiles
    HerschelBulkley.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkley_H
#define HerschelBulkley_H

#include "mixtureViscosityModel/mixtureViscosityModel.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleMultiphaseInteractingMixture;

namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HerschelBulkley Declaration
\*---------------------------------------------------------------------------*/

class HerschelBulkley
:
    public mixtureViscosityModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary HerschelBulkleyCoeffs_;

        //- HerschelBulkley viscosity coefficient
        dimensionedScalar k_;

        //- HerschelBulkley viscosity exponent
        dimensionedScalar n_;

        //- HerschelBulkley tau0
        dimensionedScalar tau0_;

        //- Minimum viscosity
        dimensionedScalar mu0_;

        //- Maximum viscosity
        dimensionedScalar muMax_;

        //- variable parameter
        bool fixedParam_;

        //- dry solids input
        bool isDrySolids_;

        //- Solid density
        dimensionedScalar rhod_;

        //- HerschelBulkley phase fraction
        const volScalarField& alpha_;



public:

    //- Runtime type information
    TypeName("HerschelBulkley");


    // Constructors

        //- Construct from components
        HerschelBulkley
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word modelName=typeName
        );


    //- Destructor
    ~HerschelBulkley()
    {}


    // Member Functions

        //- Return the strain rate mag(symm(grad(U)))
        virtual tmp<volScalarField> strainRate() const;

        //- Return the exponent
        virtual tmp<volScalarField> n() const;

        //- Return the coefficient
        virtual tmp<volScalarField> k() const;

        //- Return the yield stress
        virtual tmp<volScalarField> tau0() const;

        //- Return the concentration
        virtual tmp<volScalarField> X() const;

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        tmp<volScalarField> mu(const volScalarField& muc) const;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
