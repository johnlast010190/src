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
    Foam::homogeneousMixture

Group
    grpReactionThermophysicalMixtures

Description
    Foam::homogeneousMixture

SourceFiles
    homogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousMixture_H
#define homogeneousMixture_H

#include "mixtures/basicCombustionMixture/basicCombustionMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class homogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class homogeneousMixture
:
    public basicCombustionMixture
{
    // Private data

        static const int nSpecies_ = 1;
        static const char* specieNames_[1];

        ThermoType reactants_;
        ThermoType products_;

        mutable ThermoType mixture_;

        //- Construct as copy (not implemented)
        homogeneousMixture(const homogeneousMixture<ThermoType>&);

        //- Regress variable
        volScalarField& b_;


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, obr and phase name
        homogeneousMixture(const dictionary&, const objectRegistry&, const word&);


    //- Destructor
    virtual ~homogeneousMixture()
    {}


    // Member functions

        const ThermoType& mixture(const scalar) const;

        const ThermoType& cellMixture(const label celli) const
        {
            return mixture(b_[celli]);
        }

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture(b_.boundaryField()[patchi][facei]);
        }

        const ThermoType& cellReactants(const label) const
        {
            return reactants_;
        }

        const ThermoType& patchFaceReactants(const label, const label) const
        {
            return reactants_;
        }

        const ThermoType& cellProducts(const label) const
        {
            return products_;
        }

        const ThermoType& patchFaceProducts(const label, const label) const
        {
            return products_;
        }

        //- Read dictionary
        void read(const objectRegistry&, const dictionary&);

        //- Return thermo based on index
        const ThermoType& getLocalThermo(const label speciei) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "mixtures/homogeneousMixture/homogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
