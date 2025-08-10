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
    (c) 2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "searchableSurfaceModifier.H"
#include "triSurface/triSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceModifier, 0);
    defineRunTimeSelectionTable(searchableSurfaceModifier, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceModifier::searchableSurfaceModifier
(
    const searchableSurfaces& geometry,
    const dictionary& dict
)
:
    geometry_(geometry),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceModifier::~searchableSurfaceModifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::searchableSurfaceModifier>
Foam::searchableSurfaceModifier::New
(
    const word& type,
    const searchableSurfaces& geometry,
    const dictionary& dict
)
{
    const auto ctor = ctorTableLookup("searchableSurfaceModifier", dictionaryConstructorTable_(), type);
    return autoPtr<searchableSurfaceModifier>(ctor(geometry, dict));
}


// ************************************************************************* //
