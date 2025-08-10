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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMeshTopoChangers/none/fvMeshTopoChangersNone.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMeshTopoChanger> Foam::fvMeshTopoChanger::New
(
    fvMesh& mesh
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

    if (dictHeader.typeHeaderOk<IOdictionary>())
    {
        IOdictionary dict(dictHeader);

        // Check if the old keyword 'dynamicFvMesh' is present
        if (dict.found("dynamicFvMesh"))
        {
            FatalErrorInFunction
                << "Old keyword 'dynamicFvMesh' found in the dynamicMeshDict "
                << "dictionary. Please, make sure that the new dynamic mesh "
                << "structure is used and remove that keyword!"
                << exit(FatalError);
        }

        if (dict.found("topoChanger"))
        {
            const dictionary& topoChangerDict = dict.subDict("topoChanger");

            const word fvMeshTopoChangerTypeName
            (
                topoChangerDict.lookup("type")
            );

            Info<< "Selecting fvMeshTopoChanger "
                << fvMeshTopoChangerTypeName << endl;

            const_cast<Time&>(mesh.time()).libs().open
            (
                topoChangerDict,
                "libs",
                fvMeshConstructorTable_()
            );

            auto ctor = ctorTableLookup("fvMeshTopoChanger", fvMeshConstructorTable_(), fvMeshTopoChangerTypeName);
            return autoPtr<fvMeshTopoChanger>(ctor(mesh));
        }
    }

    return autoPtr<fvMeshTopoChanger>(new fvMeshTopoChangers::none(mesh));
}


// ************************************************************************* //
