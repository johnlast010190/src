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
    (c) 2021 Esi Ltd.

Class
    Foam::EulerImplicit

Description
    An Euler implicit solver for chemistry

SourceFiles
    EulerImplicit.C

\*---------------------------------------------------------------------------*/

#ifndef EulerImplicit_H
#define EulerImplicit_H

#include "chemistrySolver/chemistrySolver/chemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class EulerImplicit Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class EulerImplicit
:
    public chemistrySolver<ChemistryModel>
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;


        // Model constants

            //- Chemistry timescale
            scalar cTauChem_;

            //- Equilibrium rate limiter flag (on/off)
            Switch eqRateLimiter_;

        // Solver data
        mutable scalarField cTp_;

        //- Object registry
        const objectRegistry& obr_;

        //- Reference to the temperature field
        volScalarField& T_;

        //- Reference to pressure field
        volScalarField& p_;

        //- Reference to the energy field
        volScalarField& he_;

        //- Coefficients
        mutable scalarField mixCoeffs_;


public:

    //- Runtime type information
    TypeName("EulerImplicit");


    // Constructors

        //- Construct from mesh and phase name
        EulerImplicit(const objectRegistry& obr, const word& phaseName);


    //- Destructor
    virtual ~EulerImplicit();


    // Member Functions

    void updateRRInReactionI
    (
        const label index,
        const scalar pr,
        const scalar pf,
        const scalar corr,
        const label lRef,
        const label rRef,
        const scalar p,
        const scalar T,
        simpleMatrix<scalar>& RR
    ) const;

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT,
            label celli
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistrySolver/EulerImplicit/EulerImplicit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
