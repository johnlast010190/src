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
    (c) 2015 OpenCFD Ltd.
    (c) 2015-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "decompositionConstraints/decompositionConstraint/decompositionConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionConstraint, 1);
    defineRunTimeSelectionTable(decompositionConstraint, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraint::decompositionConstraint
(
    const dictionary& constraintsDict,
    const word& type
)
:
    coeffDict_(constraintsDict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionConstraint>
Foam::decompositionConstraint::New
(
    const dictionary& dict,
    const word& modelType
)
{
    Info<< "Selecting decompositionConstraint " << modelType << endl;

    const auto ctor = ctorTableLookup("decompositionConstraint type", dictionaryConstructorTable_(), modelType);
    return autoPtr<decompositionConstraint>
    (
        ctor(dict, modelType)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::decompositionConstraint::~decompositionConstraint()
{}


// ************************************************************************* //
