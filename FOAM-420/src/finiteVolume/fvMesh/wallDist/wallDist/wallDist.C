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
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "meshes/polyMesh/polyPatches/derived/baseWall/wall.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallDist, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelHashSet Foam::wallDist::findPatchIDs
(
    const fvMesh& mesh
) const
{
    labelHashSet patchIDs(mesh.boundaryMesh().findPatchIDs<wall>());

    labelHashSet patchIDsMod(mesh.boundaryMesh().size());

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        label patchi =  iter.key();
        if (mesh.boundary()[patchi].isActive())
        {
            patchIDsMod.insert(patchi);
        }
    }

    return patchIDsMod;
}


const Foam::dictionary& Foam::wallDist::setPDM
(
    const fvMesh& mesh,
    const word& patchTypeName
)
{
    dictionary fvSchemesDict = mesh.schemes().localSchemeDict();
    const word pdmType(patchTypeName+"Dist");
    if (!mesh.schemes().found(pdmType))
    {
        WarningInFunction
            << " Cannot find patch displacement method in fvSchemes."
            << " Setting to meshWave" << endl;

        fvSchemesDict.add(word(pdmType), dictionary(), false);

        fvSchemesDict.subDict(pdmType).add
        (
            word("method"),
            word("meshWave"),
            true
        );
        fvSchemes& fvs = const_cast<Foam::fvSchemes&>(mesh.schemes());
        fvs.setLocalSchemeDict(fvSchemesDict);
    }

    return mesh.schemes().subDict(pdmType);
}


void Foam::wallDist::constructn() const
{
    n_ = tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "n" & patchTypeName_,
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedVector("n" & patchTypeName_, dimless, Zero),
            patchDistMethod::patchTypes<vector>(mesh(), patchIDs_)
        )
    );

    const fvPatchList& patches = mesh().boundary();

    volVectorField::Boundary& nbf = n_.ref().boundaryFieldRef();

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();
        nbf[patchi].forceAssign(patches[patchi].nf());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDist::wallDist
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const word& patchTypeName
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, wallDist>(mesh),
    patchIDs_(findPatchIDs(mesh)),
    patchTypeName_(patchTypeName),
    dict_
    (
        static_cast<const fvSchemes&>(mesh).subOrEmptyDict
        (
            patchTypeName_ & "Dist"
        )
    ),
    pdm_
    (
        patchDistMethod::New
        (
            setPDM(mesh,patchTypeName),
            mesh,
            patchIDs_
        )
    ),
    y_
    (
        IOobject
        (
            "y" & patchTypeName_,
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y" & patchTypeName_, dimLength, SMALL),
        patchDistMethod::patchTypes<scalar>(mesh, patchIDs_)
    ),
    nRequired_(dict_.lookupOrDefault<Switch>("nRequired", false)),
    n_(volVectorField::null()),
    updateInterval_(dict_.lookupOrDefault<label>("updateInterval", 1)),
    requireUpdate_(true)
{
    if (nRequired_)
    {
        constructn();
    }

    movePoints();
}


Foam::wallDist::wallDist
(
    const fvMesh& mesh,
    const word& defaultPatchDistMethod,
    const labelHashSet& patchIDs,
    const word& patchTypeName
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, wallDist>(mesh),
    patchIDs_(patchIDs),
    patchTypeName_(patchTypeName),
    dict_
    (
        static_cast<const fvSchemes&>(mesh).subOrEmptyDict
        (
            patchTypeName_ & "Dist"
        )
    ),
    pdm_
    (
        patchDistMethod::New
        (
            setPDM(mesh,patchTypeName),
            mesh,
            patchIDs_,
            defaultPatchDistMethod
        )
    ),
    y_
    (
        IOobject
        (
            "y" & patchTypeName_,
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y" & patchTypeName_, dimLength, SMALL),
        patchDistMethod::patchTypes<scalar>(mesh, patchIDs_)
    ),
    nRequired_(dict_.lookupOrDefault<Switch>("nRequired", false)),
    n_(volVectorField::null()),
    updateInterval_(dict_.lookupOrDefault<label>("updateInterval", 1)),
    requireUpdate_(true)
{
    if (nRequired_)
    {
        constructn();
    }

    movePoints();
}


Foam::wallDist::wallDist(const fvMesh& mesh, const word& patchTypeName)
:
    wallDist
    (
        mesh,
        mesh.boundaryMesh().findPatchIDs<wall>(),
        patchTypeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDist::~wallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::wallDist::n() const
{
    if (isNull(n_()))
    {
        WarningInFunction
            << "n requested but 'nRequired' not specified in the "
            << (patchTypeName_ & "Dist") << " dictionary" << nl
            << "    Recalculating y and n fields." << endl;

        nRequired_ = true;
        constructn();
        pdm_->correct(y_, n_.ref());
    }

    return n_();
}


bool Foam::wallDist::movePoints()
{
    if
    (
        (updateInterval_ != 0)
     && ((mesh_.time().timeIndex() % updateInterval_) == 0)
    )
    {
        requireUpdate_ = true;
    }

    if (requireUpdate_ && pdm_->movePoints())
    {
        DebugInformation<< "Updating wall distance" << endl;

        requireUpdate_ = false;

        if (nRequired_)
        {
            return pdm_->correct(y_, n_.ref());
        }
        else
        {
            return pdm_->correct(y_);
        }
    }
    else
    {
        return false;
    }
}


void Foam::wallDist::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Info<< "Call to void Foam::wallDist::updateMesh(const mapPolyMesh& mpm) \n"
            << "... updating wall distance ..." << endl;
    }
    pdm_->updateMesh(mpm);

    // Force update if performing topology change
    // Note: needed?
    // - field would have been mapped, so if using updateInterval option (!= 1)
    //   live with error associated of not updating and use mapped values?
    requireUpdate_ = true;
    movePoints();
}


// ************************************************************************* //
