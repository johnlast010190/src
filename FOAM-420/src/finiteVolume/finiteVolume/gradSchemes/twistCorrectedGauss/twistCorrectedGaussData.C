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
    (c) 2011 OpenFOAM Foundation
    (c) 2016-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/twistCorrectedGauss/twistCorrectedGaussData.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(twistCorrectedGaussData, 0);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::twistCorrectedGaussData::twistCorrectedGaussData(const fvMesh& mesh)
:
    MeshObject<fvMesh, MoveableMeshObject, twistCorrectedGaussData>(mesh),
    skewCorrTensorsPtr_(nullptr),
    consistentVolsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::twistCorrectedGaussData::~twistCorrectedGaussData()
{
    deleteDemandDrivenData(skewCorrTensorsPtr_);
    deleteDemandDrivenData(consistentVolsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twistCorrectedGaussData::makeData() const
{
    if (debug)
    {
        Info<< "twistCorrectedGaussNormals::makeSkewCorrTensors() :"
            << "Constructing skew-correction tensors and consistent volumes"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    // Set local references to mesh data
    //includes boundary faces as opposed to owner()
    const labelUList& owner = mesh.faceOwner();

    const labelUList& neighbour = mesh.faceNeighbour();
    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceVectorField& Cf = mesh.Cf();

    // Set up temporary storage for edge centres and calculate
    surfaceVectorField Ce
    (
        IOobject
        (
            "Ce",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimLength
    );

    forAll(neighbour, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        Ce[facei] = w[facei]*C[own] + (1-w[facei])*C[nei];
    }

    forAll(w.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        fvsPatchVectorField& pCe = Ce.boundaryFieldRef()[patchi];
        const labelUList& faceCells = pw.patch().faceCells();

        tmp<vectorField> tpd = pw.patch().delta();
        const vectorField& pd = tpd();

        if (pw.coupled())
        {
            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];
                pCe[patchFacei] = C[own] + (1-pw[patchFacei])*pd[patchFacei];
            }
        }
        else
        {
            forAll(pw, patchFacei)
            {
                const label own = faceCells[patchFacei];

                // Boundary values are assumed specified at a point projected
                // orthogonally to the boundary from the internal point
                pCe[patchFacei] = C[own] + pd[patchFacei];
            }
        }
    }


    skewCorrTensorsPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "skewCorrTensors",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedTensor("0", dimVolume, tensor::zero)
        );
    surfaceTensorField& skewCorrTensors = *skewCorrTensorsPtr_;

    surfaceTensorField::Boundary& skewCorrTensorsbf =
        skewCorrTensors.boundaryFieldRef();

    consistentVolsPtr_ =
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "consistentVol",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimVolume, 0.0)
        );
    DimensionedField<scalar, volMesh>& consistentVol = *consistentVolsPtr_;

    const pointField& p = mesh.points();

    scalar oneThird(1.0/3.0);

    //internal face contributions
    forAll(mesh.owner(), facei)
    {
        const labelList& f = mesh.faces()[facei];
        label nPoints = f.size();

        const vector& cCf(Cf[facei]);
        const vector& cSf(Sf[facei]);

        // If the face is a triangle, use existing data
        if (nPoints == 3)
        {
            consistentVol[owner[facei]] += (cSf & (cCf-C[owner[facei]]));

            consistentVol[neighbour[facei]]
                -= (cSf & (cCf-C[neighbour[facei]]));

            skewCorrTensors[facei] = cSf * (cCf - Ce[facei]);
        }
        else
        {
            //split face into triangles (need to do manually since we must use face
            //centre decomposition
            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector subCf = oneThird*(p[f[pi]] + nextPoint + cCf);
                vector subSf = 0.5*(nextPoint - p[f[pi]])^(cCf - p[f[pi]]);

                consistentVol[owner[facei]]
                    += (subSf & (subCf-C[owner[facei]]));

                consistentVol[neighbour[facei]]
                    -= (subSf & (subCf-C[neighbour[facei]]));

                skewCorrTensors[facei] += subSf * (subCf - Ce[facei]);
            }
        }
    }

    //boundary face contributions
    //must include empty patch contributions to volume
    forAll(mesh.boundary(), pI)
    {
        //access polyPatch
        const polyPatch& pp(mesh.boundaryMesh()[pI]);

        bool indirectPatch = isA<indirectPolyPatch>(pp);

        //loop through patch faces
        const labelUList& pc(pp.faceCells());
        const vectorField::subField pCf(pp.faceCentres());
        const vectorField::subField pSf(pp.faceAreas());

        //empty switch
        bool notEmpty(mesh.boundary()[pI].size());

        forAll(pp, pfI)
        {
            const labelList& cFace = pp[pfI];
            label nPoints = cFace.size();
            label faceCell = pc[pfI];

            if (nPoints == 3) //triangular face - no decomposition needed
            {
                // Internal faces have already accounted for contributions from
                // indirect patches
                if (!indirectPatch)
                {
                    consistentVol[faceCell] += (pSf[pfI] & (pCf[pfI]-C[faceCell]));
                }

                if (notEmpty)
                {
                    vector scorr = pCf[pfI] - Ce.boundaryField()[pI][pfI];
                    skewCorrTensorsbf[pI][pfI] =
                        pSf[pfI] * scorr;
                }
            }
            else
            {
                forAll(cFace, pi)
                {
                    const point& nextPoint = p[cFace[(pi + 1) % nPoints]];
                    const point& thisPoint = p[cFace[pi]];
                    const point& faceCentre = pCf[pfI];

                    vector subCf
                        = oneThird*(thisPoint + nextPoint + faceCentre);
                    vector subSf
                        = 0.5*(nextPoint - thisPoint)^(faceCentre - thisPoint);

                    // Internal faces have already accounted for contributions from
                    // indirect patches
                    if (!indirectPatch)
                    {
                        consistentVol[faceCell] += (subSf & (subCf - C[faceCell]));
                    }

                    if (notEmpty)
                    {
                        vector scorr = subCf - Ce.boundaryField()[pI][pfI];
                        skewCorrTensorsbf[pI][pfI]
                            += subSf * scorr;
                    }

                }
            }
        }
    }

    consistentVol *= (1.0/3.0);

    fvc::applyFaceMaskTo(skewCorrTensors);

/*
    const faceList& fs = mesh.faces();
    label currPatch = -1;
    label currStartFace = 0;

    forAll(fs, facei)
    {
        while
        (
            currPatch < mesh.boundary().size() &&
            (
                (currPatch < 0 && facei >= mesh.nInternalFaces()) ||
                (
                    currPatch >= 0 &&
                    facei >= currStartFace+mesh.boundary()[currPatch].size()
                )
            )
        )
        {
            if (currPatch == mesh.boundary().size()-1)
            {
                // Must be empty patches - signal by setting invalid patch no
                currStartFace += mesh.boundary()[currPatch].size();
                currPatch++;
            }
            else
            {
                currPatch++;
                currStartFace = mesh.boundary()[currPatch].start();
            }
        }

        const labelList& f = fs[facei];
        label nPoints = f.size();

        // If the face is a triangle, use existing data
        if (nPoints == 3)
        {
            // We could just use the existing data except not accessible for
            // empty patches
            vector thisCf = (1.0/3)*(p[f[0]] + p[f[1]] + p[f[2]]);
            vector thisSf = 0.5*(p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]);

            consistentVol[owner[facei]] +=
                (1.0/3)*(thisSf & (thisCf-C[owner[facei]]));

            if (currPatch < 0) // Internal faces
            {
                consistentVol[neighbour[facei]] -=
                    (1.0/3)*(thisSf & (thisCf-C[neighbour[facei]]));
                vector scorr = thisCf - Ce[facei];
                skewCorrTensors[facei] = thisSf * scorr;
            }
            else if (currPatch < mesh.boundary().size()) // Non-empty boundaries
            {
                const label bfacei = facei-currStartFace;
                vector scorr =
                    thisCf -
                    Ce.boundaryField()[currPatch][bfacei];
                skewCorrTensors.boundaryField()[currPatch][bfacei] =
                    thisSf * scorr;
            }

        }
        else
        {
            point fCentre = p[f[0]];
            for (label pi = 1; pi < nPoints; pi++)
            {
                fCentre += p[f[pi]];
            }
            fCentre /= nPoints;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector subCf = (1.0/3)*(p[f[pi]] + nextPoint + fCentre);
                vector subSf = 0.5*(nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);

                consistentVol[owner[facei]] +=
                    (1.0/3)*(subSf & (subCf-C[owner[facei]]));

                if (currPatch < 0) // Internal faces
                {
                    consistentVol[neighbour[facei]] -=
                        (1.0/3)*(subSf & (subCf-C[neighbour[facei]]));
                    vector scorr = subCf - Ce[facei];
                    skewCorrTensors[facei] += subSf * scorr;
                }
                else if (currPatch < mesh.boundary().size()) // Non-empty boundaries
                {

                    const label bfacei = facei-currStartFace;
                    vector scorr =
                        subCf -
                        Ce.boundaryField()[currPatch][bfacei];
                    skewCorrTensors.boundaryField()[currPatch][bfacei] +=
                        subSf * scorr;
                }

            }
        }
    }
    */

    if (debug)
    {
        Info<< "twistCorrectedGaussNormals::makeCoeffs() :"
            << "Finished constructing skew Gauss data"
            << endl;
    }
}


const Foam::surfaceTensorField&
Foam::twistCorrectedGaussData::skewCorrTensors() const
{
    if (!skewCorrTensorsPtr_)
    {
        makeData();
    }

    return *skewCorrTensorsPtr_;
}

const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::twistCorrectedGaussData::consistentVols() const
{
    if (!consistentVolsPtr_)
    {
        makeData();
    }

    return *consistentVolsPtr_;
}


bool Foam::twistCorrectedGaussData::movePoints()
{
    deleteDemandDrivenData(skewCorrTensorsPtr_);
    deleteDemandDrivenData(consistentVolsPtr_);

    return true;
}


// ************************************************************************* //
