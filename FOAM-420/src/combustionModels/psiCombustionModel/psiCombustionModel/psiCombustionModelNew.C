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

#include "psiCombustionModel/psiCombustionModel/psiCombustionModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::combustionModels::psiCombustionModel>
Foam::combustionModels::psiCombustionModel::New
(
    const objectRegistry& obr,
    const word& combustionProperties,
    const word& phaseName
)
{
    const word combModelName
    (
        IOdictionary
        (
            IOobject
            (
                IOobject::groupName(combustionProperties, phaseName),
                obr.time().constant(),
                obr,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("combustionModel")
    );

    Info<< "Selecting combustion model " << combModelName << endl;

    const auto ctor = ctorTableLookup("psiCombustionModel type", dictionaryConstructorTable_(), combModelName);
    const word className = combModelName(0, combModelName.find('<'));

    return autoPtr<psiCombustionModel>
    (
        ctor(className, obr, combustionProperties, phaseName)
    );
}


// ************************************************************************* //
