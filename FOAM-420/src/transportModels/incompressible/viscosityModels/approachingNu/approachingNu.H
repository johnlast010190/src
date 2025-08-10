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
    (c) 2011 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::viscosityModels::approachingNu

Description
    A material model which is approaching the desired viscosity nu_desired,
    starting from nu_start.

    A blended algorithm is:
    if (current_time < some_fixed_time_)
    {
        scalar alpha = current_time / some_fixed_time;
        return (alpha * nu_desired_ + (1 - alpha) * nu_start_);
    }
    else
    {
        return nu_desired_;
    }

    However this one behaves better in my tests:
    if (current_time < some_fixed_time_)
    {
        scalar alpha = (nu_start * 1000 / nu_desired) ^ (1/some_fixed_time)
        return (nu_start / (alpha ^ current_time)) + nu_desired;
    }
    else
    {
        return nu_desired_;
    }


SourceFiles
    approachingNu.C

\*---------------------------------------------------------------------------*/

#ifndef approachingNu_H
#define approachingNu_H

#include "viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class approachingNu Declaration
\*---------------------------------------------------------------------------*/

class approachingNu
:
    public viscosityModel
{
    // Private data

        dictionary approachingNuCoeffs_;

        dimensionedScalar nu_desired_;
        dimensionedScalar nu_start_;
        dimensionedScalar some_fixed_time_;

        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("approachingNu");


    // Constructors

        //- construct from components
        approachingNu
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        ~approachingNu()
        {}


    // Member Functions

        //- Return the laminar viscosity
    // tmp<volScalarField> nu() const
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        // const volScalarField& nu() const { return nu_; }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// #endif
