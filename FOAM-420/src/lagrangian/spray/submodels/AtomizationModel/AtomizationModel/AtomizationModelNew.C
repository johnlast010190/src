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

\*---------------------------------------------------------------------------*/

#include "submodels/AtomizationModel/AtomizationModel/AtomizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::AtomizationModel<CloudType>>
Foam::AtomizationModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word AtomizationModelType(dict.lookup("atomizationModel"));

    Info<< "Selecting AtomizationModel " << AtomizationModelType << endl;

    const auto ctor = ctorTableLookup("AtomizationModelType type", dictionaryConstructorTable_(), AtomizationModelType);
    return autoPtr<AtomizationModel<CloudType>>(ctor(dict, owner));
}


// ************************************************************************* //
