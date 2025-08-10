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
    (c) 2022 OpenFOAM Foundation
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMeshStitchers/fvMeshStitcher/fvMeshStitcher.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshStitcher> Foam::fvMeshStitcher::New
(
    fvMesh& mesh,
    const bool changing
)
{
    IOobject dictHeader
    (
        "dynamicMeshDict",
        mesh.time().constant(),
        mesh.dbDir(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (changing && dictHeader.typeHeaderOk<IOdictionary>())
    {
        IOdictionary dict(dictHeader);

        forAllConstIter(dictionary, dict, iter)
        {
            if (!iter().isDict()) continue;

            const dictionary& changerDict = iter().dict();

            if (!changerDict.found("libs")) continue;

            const_cast<Time&>(mesh.time()).libs().open
            (
                changerDict,
                "libs",
                fvMeshConstructorTable_()
            );
        }
    }

    for (const auto& p : fvMeshConstructorTable_())
    {
        autoPtr<fvMeshStitcher> stitcherPtr(p.second(mesh));

        if
        (
            stitcherPtr->changing() ==
                (changing && dictHeader.typeHeaderOk<IOdictionary>())
        )
        {
            return stitcherPtr;
        }
    }

    FatalErrorInFunction
        << typeName << " for " << (changing ? "" : "non-")
        << "changing mesh not found " << nl << nl
        << exit(FatalError);

    return autoPtr<fvMeshStitcher>(nullptr);
}


// ************************************************************************* //
