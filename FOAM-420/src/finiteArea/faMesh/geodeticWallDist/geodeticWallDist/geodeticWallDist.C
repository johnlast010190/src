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
    (c) 2016 Esi Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "faMesh/geodeticWallDist/geodeticWallDist/geodeticWallDist.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geodeticWallDist, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::geodeticWallDist::setPDM
(
    const faMesh& mesh,
    const word& patchTypeName
)
{
    dictionary faSchemesDict = mesh.schemes().localSchemeDict();
    const word pdmType(patchTypeName+"Dist");
    if (!mesh.schemes().found(pdmType))
    {
        WarningInFunction
            << " Cannot find patch displacement method in faSchemes."
            << " Setting to Poisson"<< endl;

        faSchemesDict.add(word(pdmType), dictionary(), false);

        faSchemesDict.subDict(pdmType).add
        (
            word("method"),
            word("Poisson"),
            true
        );
        faSchemes& fas = const_cast<Foam::faSchemes&>(mesh.schemes());
        fas.setLocalSchemeDict(faSchemesDict);
    }

    return mesh.schemes().subDict(pdmType);
}


void Foam::geodeticWallDist::constructn() const
{
    n_ = tmp<areaVectorField>
    (
        new areaVectorField
        (
            IOobject
            (
                "n" & patchTypeName_,
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            mesh(),
            dimensionedVector("n" & patchTypeName_, dimless, vector::zero),
            faPatchDistMethod::patchTypes<vector>(mesh(), patchIDs_)
        )
    );

    const faPatchList& patches = mesh().boundary();

    areaVectorField::Boundary& nbf = n_.ref().boundaryFieldRef();

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();
        nbf[patchi].forceAssign(patches[patchi].edgeNormals());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geodeticWallDist::geodeticWallDist
(
    const faMesh& mesh,
    const labelHashSet& patchIDs,
    const word& patchTypeName
)
:
    MeshObject<faMesh, Foam::UpdateableMeshObject, geodeticWallDist>(mesh),
    patchIDs_(patchIDs),
    islandFaces_(labelList(0)),
    patchTypeName_(patchTypeName),
    pdm_
    (
        faPatchDistMethod::New
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
            mesh.thisDb()
        ),
        mesh,
        dimensionedScalar("y" & patchTypeName_, dimLength, SMALL),
        faPatchDistMethod::patchTypes<scalar>(mesh, patchIDs_)
    ),
    nRequired_
    (
        mesh.schemes().subDict(patchTypeName_ & "Dist").
            lookupOrDefault<Switch>("nRequired", false)
    ),
    n_(areaVectorField::null())
{
    if (nRequired_)
    {
        constructn();
    }

    movePoints();
}


Foam::geodeticWallDist::geodeticWallDist
(
    const faMesh& mesh,
    const labelHashSet& patchIDs,
    const labelList& islandFaces,
    const word& patchTypeName
)
:
    MeshObject<faMesh, Foam::UpdateableMeshObject, geodeticWallDist>(mesh),
    patchIDs_(patchIDs),
    islandFaces_(islandFaces),
    patchTypeName_(patchTypeName),
    pdm_
    (
        faPatchDistMethod::New
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
            mesh.thisDb()
        ),
        mesh,
        dimensionedScalar("y" & patchTypeName_, dimLength, SMALL),
        faPatchDistMethod::patchTypes<scalar>(mesh, patchIDs_)
    ),
    nRequired_
    (
        mesh.schemes().subDict(patchTypeName_ & "Dist")
       .lookupOrDefault<Switch>("nRequired", false)
    ),
    n_(areaVectorField::null())
{
    if (nRequired_)
    {
        constructn();
    }

    movePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::geodeticWallDist::~geodeticWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::areaVectorField& Foam::geodeticWallDist::n() const
{
    if (isNull(n_()))
    {
        WarningInFunction
            << "n requested but 'nRequired' not specified in the "
            << (patchTypeName_ & "Dist") << " dictionary" << nl
            << "    Recalculating y and n fields." << endl;

        nRequired_ = true;
        constructn();
        pdm_->correct(y_,  n_.ref());
    }

    return n_();
}


bool Foam::geodeticWallDist::movePoints()
{
    if (pdm_->movePoints())
    {
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


void Foam::geodeticWallDist::updateMesh(const mapPolyMesh& mpm)
{
    pdm_->updateMesh(mpm);
    movePoints();
}


// ************************************************************************* //
