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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "mixtures/egrMixture/egrMixture.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ThermoType>
const char* Foam::egrMixture<ThermoType>::specieNames_[3] = {"ft", "b", "egr"};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::egrMixture<ThermoType>::egrMixture
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const word& phaseName
)
:
    basicCombustionMixture
    (
        thermoDict,
        speciesTable(nSpecies_, specieNames_),
        obr,
        phaseName
    ),

    stoicRatio_("stoichiometricAirFuelMassRatio", dimless, thermoDict),

    fuel_(obr, thermoDict.subDict("fuel")),
    oxidant_(obr, thermoDict.subDict("oxidant")),
    products_(obr, thermoDict.subDict("burntProducts")),

    mixture_("mixture", fuel_),

    ft_(Y("ft")),
    b_(Y("b")),
    egr_(Y("egr"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::egrMixture<ThermoType>::mixture
(
    const scalar ft,
    const scalar b,
    const scalar egr
) const
{
    if (ft < 0.0001)
    {
        return oxidant_;
    }
    else
    {

        scalar fu = b*ft + (1.0 - b)*fres(ft, stoicRatio().value());
        scalar ox = 1 - ft - (ft - fu)*stoicRatio().value();

        fu *= (1.0 - egr);
        ox *= (1.0 - egr);

        scalar pr = 1 - fu - ox;

        mixture_ = fu*fuel_;
        mixture_ += ox*oxidant_;
        mixture_ += pr*products_;

        return mixture_;
    }
}


template<class ThermoType>
void Foam::egrMixture<ThermoType>::read(const objectRegistry& obr, const dictionary& thermoDict)
{
    thermoDict.lookup("stoichiometricAirFuelMassRatio") >> stoicRatio_;

    fuel_ = ThermoType(obr, thermoDict.subDict("fuel"));
    oxidant_ = ThermoType(obr, thermoDict.subDict("oxidant"));
    products_ = ThermoType(obr, thermoDict.subDict("burntProducts"));
}


template<class ThermoType>
const ThermoType& Foam::egrMixture<ThermoType>::getLocalThermo
(
    const label speciei
) const
{
    if (speciei == 0)
    {
        return fuel_;
    }
    else if (speciei == 1)
    {
        return oxidant_;
    }
    else if (speciei == 2)
    {
        return products_;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown specie index " << speciei << ". Valid indices are 0..2"
            << abort(FatalError);

        return fuel_;
    }
}


// ************************************************************************* //
