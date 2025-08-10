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
    (c) 2011 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingSources/faceSetBlendingSource/faceSetBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "sets/topoSetSource/topoSetSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(faceSetBlendingSource, 0);
addToRunTimeSelectionTable(blendingSource, faceSetBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faceSetBlendingSource::faceSetBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    blendingFaceSet_
    (
        mesh,
        word("faceSet") + dict.dictionaryName::name(),
        mesh.nFaces()
    )
{
    PtrList<entry> setDicts(dict.lookup("setSources"));

    forAll(setDicts, seti)
    {
        const entry& setDict = setDicts[seti];

        autoPtr<topoSetSource> faceSelector =
            topoSetSource::New(setDict.keyword(), mesh, setDict.dict());

        faceSelector->applyToSet
        (
            topoSetSource::ADD,
            blendingFaceSet_
        );
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> faceSetBlendingSource::sourceField()
{
    tmp<surfaceScalarField> faceSetField
    (
        new surfaceScalarField
        (
            IOobject
            (
                blendingFaceSet_.name()+"f",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0)
        )
    );

    labelList faces = blendingFaceSet_.toc();
    surfaceScalarField::Boundary& facesetbf = faceSetField->boundaryFieldRef();
    forAll(faces, lfI)
    {
        label facei = faces[lfI];

        if (facei < mesh_.nInternalFaces())
        {
            faceSetField->operator[](facei) = 1;
        }
        else
        {
            //handle coupled faces
            label patchi = mesh_.boundaryMesh().whichPatch(facei);

            if (mesh_.boundaryMesh()[patchi].coupled())
            {
                const polyPatch& pp = mesh_.boundaryMesh()[patchi];

                label pfI = facei - pp.start();
                facesetbf[patchi][pfI] = 1;
            }
        }

    }

    return faceSetField;
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
