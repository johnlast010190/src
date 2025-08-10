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
    (c) 2012-2016 OpenFOAM Foundation

Application
    ThermoMixture

Description

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "specie/specie.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "thermo/hConst/hConstThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/thermo/thermo.H"
#include "transport/const/constTransport.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    typedef constTransport
    <
        species::thermo
        <
            hConstThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > ThermoType;

    dictionary dict(IFstream("thermoDict")());

    ThermoType t1(dict.subDict("specie1"));
    ThermoType t2(dict.subDict("specie2"));

    Info<< "Checking Cp of mixture of hConstThermo" << endl;

    Info<< "W 1, 2, (1 + 2) = " << t1.W() << " " << t2.W() << " "
        << (t1 + t2).W() << endl;

    Info<< "Cp 1, 2, 1 + 2 = " << t1.cp(1, 1) << " " << t2.cp(1, 1) << " "
        << (t1 + t2).cp(1, 1) << endl;

    ThermoType t3(t1);
    t3 += t2;
    Info<< "Cp (1 += 2) = " << t3.cp(1, 1) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
