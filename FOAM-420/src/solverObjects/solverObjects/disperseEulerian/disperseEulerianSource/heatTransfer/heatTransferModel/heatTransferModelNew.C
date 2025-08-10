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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "heatTransferModel.H"
#include "disperseEulerian/phase/phase.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decoupledEulerian::heatTransferModel>
Foam::decoupledEulerian::heatTransferModel::New
(
    const dictionary& dict,
    const phase& phase
)
{
    word heatTransferModelType(dict.lookup("heatTransferModel"));

    Info<< "Selecting heatTransferModel for "
        << phase.name() << ": " << heatTransferModelType << endl;

    const auto ctor = ctorTableLookup("heatTransferModelType type", dictionaryConstructorTable_(), heatTransferModelType);
    return ctor(dict, phase);
}


// ************************************************************************* //
