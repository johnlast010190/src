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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 Esi Ltd.

Class
    Foam::constrainPressure

Description

SourceFiles
    constrainPressure.C

\*---------------------------------------------------------------------------*/

#ifndef constrainPressure_H
#define constrainPressure_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class NullMRF Declaration
\*---------------------------------------------------------------------------*/

class NullMRF
{
public:

    // Constructors

        //- Construct null
        NullMRF()
        {}


    // Member Functions

        //- Return the argument unchanged
        template<class Type>
        inline const Type& relative(const Type& U) const
        {
            return U;
        }

        //- Return the argument unchanged
        template<class Type>
        inline const Type& relative(const Type& U, const label patchi) const
        {
            return U;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RhoType, class RAUType, class FvOptionType>
void constrainPressure
(
    volScalarField& p,
    const RhoType& rho,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rhorAU,
    const FvOptionType& fvOptions,
    const bool includeRhoUsf = true
);

template<class RAUType>
void constrainPressure
(
    volScalarField& p,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rhorAU
);

template<class RAUType, class FvOptionType>
void constrainPressure
(
    volScalarField& p,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU,
    const FvOptionType& fvOptions
);

template<class RAUType>
void constrainPressure
(
    volScalarField& p,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/constrainPressure/constrainPressure.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
