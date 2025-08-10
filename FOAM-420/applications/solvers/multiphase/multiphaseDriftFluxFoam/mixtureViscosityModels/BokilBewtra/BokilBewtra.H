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
    (c) 2020 Esi Ltd.

Class
    Foam::mixtureViscosityModels::BokilBewtra

Description
     Viscosity correction model for a generic power-law BokilBewtra.

SourceFiles
    BokilBewtra.C

\*---------------------------------------------------------------------------*/

#ifndef BokilBewtra_H
#define BokilBewtra_H

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
                           Class BokilBewtra Declaration
\*---------------------------------------------------------------------------*/

class BokilBewtra
:
    public mixtureViscosityModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary BokilBewtraCoeffs_;

        //- BokilBewtra viscosity coefficient
        dimensionedScalar BokilBewtraViscosityCoeff_;

        //- BokilBewtra viscosity exponent
        dimensionedScalar BokilBewtraViscosityExponent_;

        //- Maximum viscosity
        dimensionedScalar muMax_;

        //- Maximum viscosity
        dimensionedScalar alphac_;

        //- BokilBewtra phase fraction
        const volScalarField& alpha_;


public:

    //- Runtime type information
    TypeName("BokilBewtra");


    // Constructors

        //- Construct from components
        BokilBewtra
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word modelName=typeName
        );


    //- Destructor
    ~BokilBewtra()
    {}


    // Member Functions

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
