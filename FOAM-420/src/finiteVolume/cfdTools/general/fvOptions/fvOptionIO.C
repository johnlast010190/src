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

#include "cfdTools/general/fvOptions/fvOption.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::option::writeHeader(Ostream& os) const
{
    os.beginBlock(name_);
}


void Foam::fv::option::writeFooter(Ostream& os) const
{
    os.endBlock();
}


void Foam::fv::option::writeData(Ostream& os) const
{
    os.writeEntry("type", type());
    os.writeEntry("active", active_);
    os << nl << indent << word(type() + "Coeffs");
    coeffs_.write(os);
}


bool Foam::fv::option::read(const dictionary& dict)
{
    dict.readIfPresent("active", active_);
    coeffs_ = dict.optionalSubDict(modelType_ + "Coeffs");

    return true;
}


// ************************************************************************* //
