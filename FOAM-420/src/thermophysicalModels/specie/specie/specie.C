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
    (c) 2016 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "specie/specie.H"
#include "global/constants/constants.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(specie, 0);
}

const Foam::fvMesh& Foam::specie::getMesh(const objectRegistry& obr) const
{
    if (isA<fvSolutionRegistry>(obr))
    {
        return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
    }
    return dynamic_cast<const fvMesh&>(obr);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie::specie(const objectRegistry& obr, const dictionary& dict)
:
    mesh_(getMesh(obr)),
    name_(dict.dictName()),
    Y_(dict.subDict("specie").lookupOrDefault("massFraction", 1.0)),
    molWeight_(readScalar(dict.subDict("specie").lookup("molWeight")))
{
    if (RR == 0.0 || Pstd == 0.0 || Tstd == 0.0)
    {
        const_cast<scalar&>(RR) = constant::physicoChemical::R.value()*1000;
        const_cast<scalar&>(Pstd) = constant::standard::Pstd.value();
        const_cast<scalar&>(Tstd) = constant::standard::Tstd.value();
        WarningInFunction
            << "Static constants not initialised!"
            << " Attmpting in-situ assignment."
            << nl << "Specie static variable values:"
            << nl << tab << "RR = " << RR
            << nl << tab << "Pstd = " << Pstd
            << nl << tab << "Tstd = " << Tstd << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::specie::write(Ostream& os) const
{
    dictionary dict("specie");
    if (Y_ != 1)
    {
        dict.add("massFraction", Y_);
    }
    dict.add("molWeight", molWeight_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie& st)
{
    st.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
