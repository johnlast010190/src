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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceReport/surfaceReport.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/nonConformalOrig/nonConformalOrigFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(surfaceReport, 0);
    addToRunTimeSelectionTable(functionObject, surfaceReport, dictionary);
}
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

void Foam::functionObjects::surfaceReport::writeFileHeader() const
{
    Log << "    Logging surface statistics to file" << endl;

    fileWriter_.writeCommented("Time");

    if (compressible_)
    {
        fileWriter_.writeDelimited("massFlux");
    }
    else
    {
        fileWriter_.writeDelimited("volumeFlux");
    }

    forAll(fields_, fI)
    {
        fileWriter_.writeDelimited(fields_[fI] + "min");
        fileWriter_.writeDelimited(fields_[fI] + "max");
        fileWriter_.writeDelimited(fields_[fI] + "mean");

        if (addDirectionMeanOutput_)
        {
            fileWriter_.writeDelimited(fields_[fI] + "directionMean");
        }

        if (homogeneity_)
        {
            fileWriter_.writeDelimited(fields_[fI] + "homogeneity");
        }
        else
        {
            fileWriter_.writeDelimited(fields_[fI] + "stDev");
        }
    }

    if (swirl_)
    {
        fileWriter_.writeDelimited("swirl");
    }
    if (tumble_)
    {
        fileWriter_.writeDelimited("tumble");
    }
    if (backFlowReport_)
    {
        fileWriter_.writeDelimited("backFlowPercentage");
    }

    fileWriter_.endLine();
}


void functionObjects::surfaceReport::reduceMaxVectorMag(vector& v)
{
    List<vector> vl(Pstream::nProcs(), vector::zero);

    vl[Pstream::myProcNo()] = v;

    Pstream::allGatherList(vl);

    v = vl[0];
    scalar vmag = mag(v);

    forAll(vl, procI)
    {
        scalar cvMag = mag(vl[procI]);
        if (cvMag > vmag)
        {
            vmag = cvMag;
            v = vl[procI];
        }
    }
}


Foam::tmp<Foam::surfaceScalarField>
functionObjects::surfaceReport::correctedFlux()
{
    return fvc::applyFaceMask(lookupObject<surfaceScalarField>(fluxName_));
}


void functionObjects::surfaceReport::totalFlux()
{
    totalFlux_ = 0;

    tmp<surfaceScalarField> tflux = correctedFlux();
    const surfaceScalarField& flux = tflux();

    forAll(flux, rfI)
    {
        if (reportFaces_.get(rfI))
        {
            totalFlux_ += mag(flux[rfI]);
        }
    }

    forAll(flux.boundaryField(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];
        const label pstart = pp.start();

        if (!isA<nonConformalOrigPolyPatch>(pp))
        {
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    totalFlux_ += mag(flux.boundaryField()[pI][pfI]);
                }
            }
        }
        else
        {
            // Add total fluxes from non-conformal patches
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    const label polyBFacei =
                        pstart + pfI - mesh_.nInternalFaces();
                    const labelUList patches =
                        mesh_.polyBFacePatches()[polyBFacei];
                    const labelUList patchFaces =
                        mesh_.polyBFacePatchFaces()[polyBFacei];

                    forAll(patches, i)
                    {
                        totalFlux_ +=
                            mag
                            (
                                flux.boundaryField()[patches[i]][patchFaces[i]]
                            );
                    }
                }
            }
        }
    }

    reduce(totalFlux_, sumOp<scalar>());
}


void functionObjects::surfaceReport::relativeFlux()
{
    relativeFlux_ = 0;

    tmp<surfaceScalarField> tflux = correctedFlux();
    const surfaceScalarField& flux = tflux();
    const surfaceVectorField& Sf = mesh_.Sf();

    forAll(flux, rfI)
    {
        if (rfI < mesh_.nInternalFaces())
        {
            if (reportFaces_.get(rfI))
            {
                relativeFlux_ +=
                    sign(Sf[rfI] & surfaceDirectionPtr_())*flux[rfI];
            }
        }
    }

    forAll(flux.boundaryField(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];
        const label pstart = pp.start();

        if (isA<processorPolyPatch>(pp))
        {
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    relativeFlux_ +=
                        sign
                        (
                            Sf.boundaryField()[pI][pfI]
                          & surfaceDirectionPtr_()
                        )
                       *flux.boundaryField()[pI][pfI];
                }
            }
        }
        else if (!isA<nonConformalOrigPolyPatch>(pp))
        {
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    relativeFlux_ += flux.boundaryField()[pI][pfI];
                }
            }
        }
        else
        {
            // Add relative fluxes from non-conformal patches
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    const label polyBFacei =
                        pstart + pfI - mesh_.nInternalFaces();
                    const labelUList patches =
                        mesh_.polyBFacePatches()[polyBFacei];
                    const labelUList patchFaces =
                        mesh_.polyBFacePatchFaces()[polyBFacei];

                    forAll(patches, i)
                    {
                        relativeFlux_ +=
                            flux.boundaryField()[patches[i]][patchFaces[i]];
                    }
                }
            }
        }
    }

    reduce(relativeFlux_, sumOp<scalar>());
}


void functionObjects::surfaceReport::backFlowArea()
{
    backFlowArea_ = 0.0;

    tmp<surfaceScalarField> tflux = correctedFlux();
    const surfaceScalarField& flux = tflux();
    const surfaceVectorField& Sf = mesh_.Sf();

    forAll(flux, rfI)
    {
        if (rfI < mesh_.nInternalFaces())
        {
            if (reportFaces_.get(rfI))
            {
                if (sign(Sf[rfI] & surfaceDirectionPtr_())*flux[rfI] < 0)
                {
                    backFlowArea_ += mesh_.magSf()[rfI];
                }
            }
        }
    }

    forAll(flux.boundaryField(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];
        const label pstart = pp.start();

        if (isA<processorPolyPatch>(pp))
        {
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    if
                    (
                        sign
                        (
                            Sf.boundaryField()[pI][pfI]
                          & surfaceDirectionPtr_()
                        )*flux.boundaryField()[pI][pfI] < 0
                    )
                    {
                        backFlowArea_ += mesh_.magSf().boundaryField()[pI][pfI];
                    }
                }
            }
        }
        else if (!isA<nonConformalOrigPolyPatch>(pp))
        {
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    if (flux.boundaryField()[pI][pfI] < 0)
                    {
                        backFlowArea_ += mesh_.magSf().boundaryField()[pI][pfI];
                    }
                }
            }
        }
        else
        {
            // Add backflow area from non-conformal patches
            forAll(flux.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    const label polyBFacei =
                        pstart + pfI - mesh_.nInternalFaces();
                    const labelUList patches =
                        mesh_.polyBFacePatches()[polyBFacei];
                    const labelUList patchFaces =
                        mesh_.polyBFacePatchFaces()[polyBFacei];

                    forAll(patches, i)
                    {
                        if (flux.boundaryField()[patches[i]][patchFaces[i]] < 0)
                        {
                            backFlowArea_ +=
                                mesh_.magSf().boundaryField()
                                [patches[i]][patchFaces[i]];
                        }
                    }
                }
            }
        }
    }

    reduce(backFlowArea_, sumOp<scalar>());
}


template<class GeoField>
void functionObjects::surfaceReport::calcNonScalarData
(
    const GeoField& f,
    FixedList<scalar, 5>& data
)
{
    tmp<surfaceScalarField> tflux = correctedFlux();
    const surfaceScalarField& flux = tflux();

    data[0] = GREAT;
    data[1] = -GREAT;
    data[2] = 0;
    data[3] = 0;
    data[4] = -1;    // Directional flux mean not defined

    typedef GeometricField
    <
        typename GeoField::value_type,Foam::fvsPatchField, Foam::surfaceMesh
    > SurfaceFieldType;

    autoPtr<SurfaceFieldType> ff;
    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        ff.reset
        (
            new SurfaceFieldType(fvc::applyFaceMask(fvc::interpolate(f)))
        );
    }

    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        forAll(flux, nI)
        {
            if (reportFaces_.get(nI))
            {
                data[0] = min(data[0], mag(ff()[nI]));
                data[1] = max(data[1], mag(ff()[nI]));

                if (fluxWeighting_)
                {
                    data[2] += mag(flux[nI]*ff()[nI]);
                }
                else
                {
                    data[2] += mesh_.magSf()[nI]*mag(ff()[nI]);
                }
            }
        }
    }

    forAll(mesh_.boundaryMesh(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];
        const label pstart = pp.start();

        if (!isA<nonConformalOrigPolyPatch>(pp))
        {
            const Field<typename GeoField::value_type>& pfv
            (
                haveIndirectPatches_
              ? Field<typename GeoField::value_type>(ff().boundaryField()[pI])
              : Field<typename GeoField::value_type>(f.boundaryField()[pI])
            );

            forAll(pfv, pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    data[0] = min(data[0], mag(pfv[pfI]));
                    data[1] = max(data[1], mag(pfv[pfI]));

                    if (fluxWeighting_)
                    {
                        data[2] +=
                            mag(flux.boundaryField()[pI][pfI]*pfv[pfI]);
                    }
                    else
                    {
                        data[2] +=
                            mesh_.boundary()[pI].magSf()[pfI]*mag(pfv[pfI]);
                    }
                }
            }
        }
        else
        {
            forAll(f.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    const label polyBFacei =
                        pstart + pfI - mesh_.nInternalFaces();
                    const labelUList patches =
                        mesh_.polyBFacePatches()[polyBFacei];
                    const labelUList patchFaces =
                        mesh_.polyBFacePatchFaces()[polyBFacei];

                    forAll(patches, i)
                    {
                        data[0] = min
                        (
                            data[0],
                            mag(f.boundaryField()[patches[i]][patchFaces[i]])
                        );

                        data[1] = max
                        (
                            data[1],
                            mag(f.boundaryField()[patches[i]][patchFaces[i]])
                        );

                        if (fluxWeighting_)
                        {
                            data[2] +=
                                mag
                                (
                                    flux.boundaryField()
                                    [patches[i]][patchFaces[i]]
                                   *f.boundaryField()[patches[i]][patchFaces[i]]
                                );
                        }
                        else
                        {
                            data[2] +=
                                mesh_.boundary()[patches[i]]
                               .magSf()[patchFaces[i]]*mag(f.boundaryField()
                                [patches[i]][patchFaces[i]]);
                        }
                    }
                }
            }
        }
    }

    reduce(data[0], minOp<scalar>());
    reduce(data[1], maxOp<scalar>());
    reduce(data[2], sumOp<scalar>());

    if (fluxWeighting_)
    {
        data[2] /= (totalFlux_ + SMALL);
    }
    else
    {
        data[2] /= (totalSurface_ + SMALL);
    }

    // Calculate standard deviation
    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        forAll(flux, nI)
        {
            if (reportFaces_.get(nI))
            {
                if (fluxWeighting_)
                {
                    if (homogeneity_)
                    {
                        data[3] +=
                            mag
                            (
                                flux[nI])*(1 - 0.5*(mag(1 - (mag(ff()[nI])
                               /(data[2] + SMALL))))
                            );
                    }
                    else
                    {
                        data[3] += mag(flux[nI])*sqr(mag(ff()[nI]) - data[2]);
                    }
                }
                else
                {
                    if (homogeneity_)
                    {
                        data[3] +=
                            mesh_.magSf()[nI]
                           *(1 - 0.5*(mag(1 - (mag(ff()[nI])
                                /(data[2] + SMALL)))));
                    }
                    else
                    {
                        data[3] +=
                            mesh_.magSf()[nI]*sqr(mag(ff()[nI]) - data[2]);
                    }
                }
            }
        }
    }

    forAll(mesh_.boundaryMesh(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];
        const label pstart = pp.start();

        if (!isA<nonConformalOrigPolyPatch>(pp))
        {
            const Field<typename GeoField::value_type>& pfv
            (
                haveIndirectPatches_
              ? Field<typename GeoField::value_type>(ff().boundaryField()[pI])
              : Field<typename GeoField::value_type>(f.boundaryField()[pI])
            );

            forAll(pfv, pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    if (fluxWeighting_)
                    {
                        if (homogeneity_)
                        {
                            data[3] +=
                                mag(flux.boundaryField()[pI][pfI])
                               *(1 - 0.5*
                                    (mag(1 - (mag(pfv[pfI])
                                   /(data[2] + SMALL)))));
                        }
                        else
                        {
                            data[3] +=
                                mag(flux.boundaryField()[pI][pfI])
                               *sqr(mag(pfv[pfI]) - data[2]);
                        }
                    }
                    else
                    {
                        if (homogeneity_)
                        {
                            data[3] +=
                                mesh_.boundary()[pI].magSf()[pfI]
                               *(1 - 0.5*(mag(1 - (mag(pfv[pfI])
                                   /(data[2] + SMALL)))));
                        }
                        else
                        {
                            data[3] +=
                                mesh_.boundary()[pI].magSf()[pfI]
                               *sqr(mag(pfv[pfI]) - data[2]);
                        }
                    }
                }
            }
        }
        else
        {
            forAll(f.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    const label polyBFacei =
                        pstart + pfI - mesh_.nInternalFaces();
                    const labelUList patches =
                        mesh_.polyBFacePatches()[polyBFacei];
                    const labelUList patchFaces =
                        mesh_.polyBFacePatchFaces()[polyBFacei];

                    forAll(patches, i)
                    {
                        if (fluxWeighting_)
                        {
                            if (homogeneity_)
                            {
                                data[3] +=
                                    mag
                                    (
                                        flux.boundaryField()
                                        [patches[i]][patchFaces[i]]
                                    )
                                   *(1 - 0.5*(mag(1 - (mag(f.boundaryField()
                                        [patches[i]][patchFaces[i]])
                                       /(data[2] + SMALL))))
                                    );
                            }
                            else
                            {
                                data[3] +=
                                    mag
                                    (
                                        flux.boundaryField()
                                        [patches[i]][patchFaces[i]]
                                    )
                                   *sqr(mag(f.boundaryField()
                                        [patches[i]][patchFaces[i]]) - data[2]);
                            }
                        }
                        else
                        {
                            if (homogeneity_)
                            {
                                data[3] +=
                                    mesh_.boundary()[patches[i]]
                                   .magSf()[patchFaces[i]]
                                   *(1 - 0.5*(mag(1 - (mag(f.boundaryField()
                                    [patches[i]][patchFaces[i]])
                                   /(data[2] + SMALL)))));
                            }
                            else
                            {
                                data[3] +=
                                    mesh_.boundary()[patches[i]]
                                   .magSf()[patchFaces[i]]
                                   *sqr(mag(f.boundaryField()
                                    [patches[i]][patchFaces[i]]) - data[2]);
                            }
                        }
                    }
                }
            }
        }
    }

    reduce(data[3], sumOp<scalar>());

    if (fluxWeighting_)
    {
        data[3] /= (totalFlux_ + SMALL);
    }
    else
    {
        data[3] /= (totalSurface_ + SMALL);
    }

    if (!homogeneity_)
    {
        data[3] = sqrt(mag(data[3]));
    }
}


void functionObjects::surfaceReport::calcScalarData
(
    const volScalarField& f,
    FixedList<scalar, 5>& data
)
{
    tmp<surfaceScalarField> tflux = correctedFlux();
    const surfaceScalarField& flux = tflux();

    data[0] = GREAT;     // min
    data[1] = -GREAT;    // max
    data[2] = 0;         // mean
    data[3] = 0;         // stdDev / homogeneity
    data[4] = 0;         // directionMean

    autoPtr<surfaceScalarField> ff;

    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        ff.reset
        (
            new surfaceScalarField(fvc::applyFaceMask(fvc::interpolate(f)))
        );

        forAll(flux, nI)
        {
            if (reportFaces_.get(nI))
            {
                data[0] = min(data[0], ff()[nI]);
                data[1] = max(data[1], ff()[nI]);

                if (fluxWeighting_)
                {
                    if (directionAware_)
                    {
                        data[2] += flux[nI]*ff()[nI];
                    }
                    else
                    {
                        data[2] += mag(flux[nI])*ff()[nI];
                    }

                    if (addDirectionMeanOutput_)
                    {
                        data[4] += flux[nI]*ff()[nI];
                    }
                }
                else
                {
                    data[2] += mesh_.magSf()[nI]*ff()[nI];
                }
            }
        }
    }

    forAll(mesh_.boundaryMesh(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];
        const label pstart = pp.start();

        if (!isA<nonConformalOrigPolyPatch>(pp))
        {
            const scalarField& pfv
            (
                haveIndirectPatches_
              ? scalarField(ff().boundaryField()[pI])
              : scalarField(f.boundaryField()[pI])
            );

            forAll(pfv, pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    data[0] = min(data[0], pfv[pfI]);
                    data[1] = max(data[1], pfv[pfI]);

                    if (fluxWeighting_)
                    {
                        if (directionAware_)
                        {
                            data[2] += flux.boundaryField()[pI][pfI]*pfv[pfI];
                        }
                        else
                        {
                            data[2] +=
                                mag(flux.boundaryField()[pI][pfI])*pfv[pfI];
                        }

                        if (addDirectionMeanOutput_)
                        {
                            data[4] += flux.boundaryField()[pI][pfI]*pfv[pfI];
                        }
                    }
                    else
                    {
                        data[2] += mesh_.boundary()[pI].magSf()[pfI]*pfv[pfI];
                    }
                }
            }
        }
        else
        {
            forAll(f.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    const label polyBFacei =
                        pstart + pfI - mesh_.nInternalFaces();
                    const labelUList patches =
                        mesh_.polyBFacePatches()[polyBFacei];
                    const labelUList patchFaces =
                        mesh_.polyBFacePatchFaces()[polyBFacei];

                    forAll(patches, i)
                    {
                        data[0] = min
                        (
                            data[0],
                            f.boundaryField()[patches[i]][patchFaces[i]]
                        );

                        data[1] = max
                        (
                            data[1],
                            f.boundaryField()[patches[i]][patchFaces[i]]
                        );

                        if (fluxWeighting_)
                        {
                            if (directionAware_)
                            {
                                data[2] +=
                                    flux.boundaryField()
                                    [patches[i]][patchFaces[i]]
                                   *f.boundaryField()
                                    [patches[i]][patchFaces[i]];
                            }
                            else
                            {
                                data[2] +=
                                    mag
                                    (
                                        flux.boundaryField()
                                        [patches[i]][patchFaces[i]]
                                    )
                                   *f.boundaryField()
                                    [patches[i]][patchFaces[i]];
                            }

                            if (addDirectionMeanOutput_)
                            {
                                data[4] +=
                                    flux.boundaryField()
                                    [patches[i]][patchFaces[i]]
                                   *f.boundaryField()
                                    [patches[i]][patchFaces[i]];
                            }
                        }
                        else
                        {
                            data[2] +=
                                mesh_.boundary()[patches[i]]
                               .magSf()[patchFaces[i]]
                               *f.boundaryField()[patches[i]][patchFaces[i]];
                        }
                    }
                }
            }
        }
    }

    reduce(data[0], minOp<scalar>());
    reduce(data[1], maxOp<scalar>());
    reduce(data[2], sumOp<scalar>());

    if (addDirectionMeanOutput_)
    {
        reduce(data[4], sumOp<scalar>());
    }

    if (fluxWeighting_)
    {
        if (directionAware_)
        {
            data[2] /= (relativeFlux_ + SMALL);
        }
        else
        {
            data[2] /= (totalFlux_ + SMALL);
        }

        if (addDirectionMeanOutput_)
        {
            data[4] /= (relativeFlux_ + SMALL);
        }
    }
    else
    {
        data[2] /= (totalSurface_ + SMALL);
    }

    // Calculate standard deviation
    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        forAll(flux, nI)
        {
            if (reportFaces_.get(nI))
            {
                if (fluxWeighting_)
                {
                    if (homogeneity_)
                    {
                        data[3] +=
                            mag(flux[nI])
                           *(1 - 0.5*(mag(1 - (ff()[nI]/(data[2] + SMALL)))));
                    }
                    else
                    {
                        data[3] += mag(flux[nI])*sqr(ff()[nI] - data[2]);
                    }
                }
                else
                {
                    if (homogeneity_)
                    {
                        data[3] +=
                            mesh_.magSf()[nI]
                           *(1 - 0.5*(mag(1 - (ff()[nI]/(data[2] + SMALL)))));
                    }
                    else
                    {
                        data[3] += mesh_.magSf()[nI]*sqr(ff()[nI] - data[2]);
                    }
                }
            }
        }
    }

    forAll(mesh_.boundaryMesh(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];
        const label pstart = pp.start();

        if (!isA<nonConformalOrigPolyPatch>(pp))
        {
            const scalarField& pfv
            (
                haveIndirectPatches_
              ? scalarField(ff().boundaryField()[pI])
              : scalarField(f.boundaryField()[pI])
            );

            forAll(pfv, pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    if (fluxWeighting_)
                    {
                        if (homogeneity_)
                        {
                            data[3] +=
                                mag(pfv[pfI])
                               *(1 - 0.5*(mag(1-(pfv[pfI]/(data[2] + SMALL)))));
                        }
                        else
                        {
                            data[3] +=
                                mag(flux.boundaryField()[pI][pfI])
                               *sqr(pfv[pfI] - data[2]);
                        }
                    }
                    else
                    {
                        if (homogeneity_)
                        {
                            data[3] +=
                                mesh_.boundary()[pI].magSf()[pfI]
                               *(1 - 0.5*(mag(1-(pfv[pfI]/(data[2] + SMALL)))));
                        }
                        else
                        {
                            data[3] +=
                                mesh_.boundary()[pI].magSf()[pfI]
                               *sqr(pfv[pfI] - data[2]);
                        }
                    }
                }
            }
        }
        else
        {
            forAll(f.boundaryField()[pI], pfI)
            {
                if (reportFaces_.get(pstart + pfI))
                {
                    const label polyBFacei =
                        pstart + pfI - mesh_.nInternalFaces();
                    const labelUList patches =
                        mesh_.polyBFacePatches()[polyBFacei];
                    const labelUList patchFaces =
                        mesh_.polyBFacePatchFaces()[polyBFacei];

                    forAll(patches, i)
                    {
                        if (fluxWeighting_)
                        {
                            if (homogeneity_)
                            {
                                data[3] +=
                                    mag
                                    (
                                        f.boundaryField()
                                        [patches[i]][patchFaces[i]]
                                    )
                                   *(
                                        1 - 0.5*(mag(1 - (f.boundaryField()
                                        [patches[i]][patchFaces[i]]
                                        /(data[2] + SMALL))))
                                    );
                            }
                            else
                            {
                                data[3] +=
                                    mag
                                    (
                                        flux.boundaryField()
                                        [patches[i]][patchFaces[i]]
                                    )
                                   *sqr
                                    (
                                        f.boundaryField()
                                        [patches[i]][patchFaces[i]] - data[2]
                                    );
                            }
                        }
                        else
                        {
                            if (homogeneity_)
                            {
                                data[3] +=
                                    mesh_.boundary()[patches[i]]
                                   .magSf()[patchFaces[i]]
                                   *(
                                        1 - 0.5*(mag(1 - (f.boundaryField()
                                        [patches[i]][patchFaces[i]]
                                       /(data[2] + SMALL))))
                                    );
                            }
                            else
                            {
                                data[3] +=
                                    mesh_.boundary()[patches[i]]
                                   .magSf()[patchFaces[i]]
                                   *sqr
                                    (
                                        f.boundaryField()
                                        [patches[i]][patchFaces[i]] - data[2]
                                    );
                            }
                        }
                    }
                }
            }
        }
    }

    reduce(data[3], sumOp<scalar>());

    if (fluxWeighting_)
    {
        data[3] /= (totalFlux_ + SMALL);
    }
    else
    {
        data[3] /= (totalSurface_ + SMALL);
    }

    if (!homogeneity_)
    {
        data[3] = sqrt(data[3]);
    }
}


void functionObjects::surfaceReport::calcSwirlData(const volVectorField& f)
{
    tmp<surfaceScalarField> tflux = correctedFlux();
    const surfaceScalarField& flux = tflux();

    swirlData_ = 0;

    autoPtr<surfaceVectorField> ff;

    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        ff.reset
        (
            new surfaceVectorField(fvc::applyFaceMask(fvc::interpolate(f)))
        );

        forAll(flux, nI)
        {
            if (reportFaces_.get(nI))
            {
                point fc = mesh_.faceCentres()[nI];
                point hitPoint =
                    swirlCentre_
                  + ((fc - swirlCentre_) & swirlUnitNorm_)*swirlUnitNorm_;
                point hitVec = fc - hitPoint;
                scalar dist = mag(hitVec);

                point tangentialDir = swirlUnitNorm_ ^ (hitVec/(dist + SMALL));

                swirlData_ += flux[nI]*(ff()[nI] & tangentialDir)*dist;
            }
        }
    }

    forAll(mesh_.boundaryMesh(), pI)
    {
        label pstart = mesh_.boundaryMesh()[pI].start();

        const vectorField& pfv
        (
            haveIndirectPatches_
          ? vectorField(ff().boundaryField()[pI])
          : vectorField(f.boundaryField()[pI])
        );

        forAll(pfv, pfI)
        {
            if (reportFaces_.get(pstart + pfI))
            {
                point fc = mesh_.faceCentres()[pstart + pfI];
                point hitPoint =
                    swirlCentre_
                  + ((fc - swirlCentre_) & swirlUnitNorm_)*swirlUnitNorm_;
                point hitVec = fc - hitPoint;
                scalar dist = mag(hitVec);

                point tangentialDir = swirlUnitNorm_ ^ (hitVec/(dist + SMALL));

                swirlData_ +=
                    flux.boundaryField()[pI][pfI]
                   *(pfv[pfI] & tangentialDir)*dist;
            }
        }
    }

    reduce(swirlData_, sumOp<scalar>());
}


Foam::tmp<Foam::volScalarField> functionObjects::surfaceReport::getRho() const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(lookupObject<volScalarField>(rhoName_));
    }
}


void functionObjects::surfaceReport::calcTumble(const volVectorField& U)
{
    const volScalarField rho(getRho());

    tumbleValue_ = 0;

    scalar omegaTumble = 0;
    scalar otNum = 0;    // Omega tumble numerator
    scalar otDen = 0;    // Omega tumble denominator
    scalar omegaEngine = 0;

    scalar planeArea = 0;
    scalar massFlow = 0;
    scalar UzAveraged = 0;
    scalar rhoAveraged = 0;

    autoPtr<surfaceVectorField> Uf;
    autoPtr<surfaceScalarField> rhof;

    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        Uf.reset
        (
            new surfaceVectorField(fvc::applyFaceMask(fvc::interpolate(U)))
        );

        rhof.reset
        (
            new surfaceScalarField(fvc::applyFaceMask(fvc::interpolate(rho)))
        );
    }

    if (reportContainsInternalFaces_ || haveIndirectPatches_)
    {
        // Calculate average quantities first
        forAll(Uf(), nI)
        {
            if (reportFaces_.get(nI))
            {
                scalar Uz = Uf()[nI] & tumbleDir_;
                // if (!compressible_) Uz /= rhof()[nI];
                const scalar& area = mesh_.magSf()[nI];
                scalar flux = rhof()[nI]*Uz*area;

                planeArea += area;
                massFlow += flux;

                if (massAveraged_)
                {
                    UzAveraged += Uz*flux;
                    rhoAveraged += rhof()[nI]*flux;
                }
                else
                {
                    UzAveraged += Uz*area;
                    rhoAveraged += rhof()[nI]*area;
                }
            }
        }

        reduce(planeArea, sumOp<scalar>());
        reduce(massFlow, sumOp<scalar>());

        reduce(UzAveraged, sumOp<scalar>());
        reduce(rhoAveraged, sumOp<scalar>());

        if (massAveraged_)
        {
            UzAveraged /= massFlow;
            rhoAveraged /= massFlow;
        }
        else
        {
            UzAveraged /= planeArea;
            rhoAveraged /= planeArea;
        }

        // Calculate the tumble quantities
        forAll(Uf(), nI)
        {
            if (reportFaces_.get(nI))
            {
                const point& fc = mesh_.faceCentres()[nI];
                scalar xCoor = fc.x() - tumbleOrig_.x();
                const scalar& area = mesh_.magSf()[nI];
                // const scalar& area = mesh_.faceAreas()[nI];
                scalar Uz = Uf()[nI] & tumbleDir_;
                // if (!compressible_) Uz /= rhof()[nI];

                otNum += rhof()[nI]*(Uz - UzAveraged)*xCoor*area;
                otDen += rhof()[nI]*xCoor*xCoor*area;
            }
        }

        reduce(otNum, sumOp<scalar>());
        reduce(otDen, sumOp<scalar>());
    }
    else
    {
        // Calculate average quantities first
        forAll(U.boundaryField(), pI)
        {
            label pstart = mesh_.boundaryMesh()[pI].start();

            const vectorField& Upfv
            (
                haveIndirectPatches_
              ? vectorField(Uf().boundaryField()[pI])
              : vectorField(U.boundaryField()[pI])
            );

            const scalarField& rhopfv
            (
                haveIndirectPatches_
              ? scalarField(rhof().boundaryField()[pI])
              : scalarField(rho.boundaryField()[pI])
            );

            forAll(Upfv, pfI)
            {
                label nI = pstart + pfI;
                if (reportFaces_.get(nI))
                {
                    scalar Uz = Upfv[pfI] & tumbleDir_;
                    // if (!compressible_) Uz /= rhof[nI];
                    const scalar& area = mesh_.magSf()[nI];
                    scalar flux = rhopfv[pfI]*Uz*area;

                    planeArea += area;
                    massFlow += flux;

                    if (massAveraged_)
                    {
                        UzAveraged += Uz*flux;
                        rhoAveraged += rhopfv[pfI]*flux;
                    }
                    else
                    {
                        UzAveraged += Uz*area;
                        rhoAveraged += rhopfv[pfI]*area;
                    }
                }
            }
        }

        reduce(planeArea, sumOp<scalar>());
        reduce(massFlow, sumOp<scalar>());

        reduce(UzAveraged, sumOp<scalar>());
        reduce(rhoAveraged, sumOp<scalar>());

        if (massAveraged_)
        {
            UzAveraged /= massFlow;
            rhoAveraged /= massFlow;
        }
        else
        {
            UzAveraged /= planeArea;
            rhoAveraged /= planeArea;
        }

        // Calculate the tumble quantities
        forAll(mesh_.boundaryMesh(), pI)
        {
            const vectorField& Upfv
            (
                haveIndirectPatches_
              ? vectorField(Uf().boundaryField()[pI])
              : vectorField(U.boundaryField()[pI])
            );

            const scalarField& rhopfv
            (
                haveIndirectPatches_
              ? scalarField(rhof().boundaryField()[pI])
              : scalarField(rho.boundaryField()[pI])
            );

            label pstart = mesh_.boundaryMesh()[pI].start();

            forAll(Upfv, pfI)
            {
                label nI = pstart + pfI;
                if (reportFaces_.get(nI))
                {
                    const point& fc = mesh_.faceCentres()[nI];
                    scalar xCoor = fc.x() - tumbleOrig_.x();
                    const scalar& area = mesh_.magSf()[nI];
                    // const scalar& area = mesh_.faceAreas()[nI];
                    scalar Uz = Upfv[pfI] & tumbleDir_;
                    // if (!compressible_) Uz /= rhof[nI];

                    otNum += rhopfv[pfI]*(Uz - UzAveraged)*xCoor*area;
                    otDen += rhopfv[pfI]*xCoor*xCoor*area;
                }
            }
        }

        reduce(otNum, sumOp<scalar>());
        reduce(otDen, sumOp<scalar>());
    }

    omegaTumble = otNum/otDen;
    omegaEngine =
        (constant::mathematical::pi*massFlow)
       /(Ls_*planeArea*rhoAveraged);

    tumbleValue_ = omegaTumble/omegaEngine;
}


void functionObjects::surfaceReport::markFaceZone(const dictionary& sD)
{
    word fzName(sD.lookup("name"));
    label faceZoneID = mesh_.faceZones().findZoneID(fzName);

    if (faceZoneID != -1)
    {
        const faceZone& fz = mesh_.faceZones()[faceZoneID];
        const polyPatchList& patches = mesh_.boundaryMesh();
        forAll(fz, zfI)
        {
            label gfI = fz[zfI];

            if (gfI < mesh_.nInternalFaces())
            {
                reportContainsInternalFaces_ = true;
                reportFaces_.set(gfI, 1);
            }
            else
            {
                label patchIndex = mesh_.boundaryMesh().whichPatch(gfI);

                if (isA<processorPolyPatch>(patches[patchIndex]))
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(patches[patchIndex]);

                    // Only add lower numbered proc face to avoid duplication
                    if (procPatch.myProcNo() < procPatch.neighbProcNo())
                    {
                        reportFaces_.set(gfI, 1);
                    }
                }
                else
                {
                    reportFaces_.set(gfI, 1);
                }
            }
        }
    }
    else
    {
        WarningInFunction
            << "Face zone " << fzName << " could not be found"
            << " and will not be added to " << name()
            << " surface report." << endl;
    }
}


void functionObjects::surfaceReport::markCellZonePair(const dictionary& sD)
{
    word zone1name(sD.lookup("zone1"));
    word zone2name(sD.lookup("zone2"));

    label zone1ID = mesh_.cellZones().findZoneID(zone1name);
    label zone2ID = mesh_.cellZones().findZoneID(zone2name);

    if (zone1ID != -1 && zone2ID != -1)
    {
        reportContainsInternalFaces_ = true;

        const cellZone& zone1(mesh_.cellZones()[zone1ID]);
        const cellZone& zone2(mesh_.cellZones()[zone2ID]);

        const labelList& neis = mesh_.neighbour();
        const labelList& owns = mesh_.faceOwner();

        forAll(neis, nI)
        {
            label neiCellI = neis[nI];
            label ownCellI = owns[nI];

            if
            (
                zone1.whichCell(ownCellI) != -1
             && zone2.whichCell(neiCellI) != -1
            )
            {
                reportFaces_.set(nI, 1);
            }
            else if
            (
                zone1.whichCell(neiCellI) != -1
             && zone2.whichCell(ownCellI) != -1
            )
            {
                reportFaces_.set(nI, 1);
            }
        }

        // Now check all processor faces

        // Calculate coupled zoneID
        labelList bfZoneID(mesh_.nFaces() - mesh_.nInternalFaces());

        for
        (
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            bfZoneID[faceI - mesh_.nInternalFaces()] =
                mesh_.cellZones().whichZone(owns[faceI]);
        }

        labelList coupleZoneID = bfZoneID;

        syncTools::swapBoundaryFaceList(mesh_, coupleZoneID);

        const polyPatchList& patches = mesh_.boundaryMesh();

        forAll(patches, patchI)
        {
            if (isA<processorPolyPatch>(patches[patchI]))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Only add lower numbered proc face to avoid duplication
                if (procPatch.myProcNo() < procPatch.neighbProcNo())
                {
                    forAll(procPatch, pfI)
                    {
                        label faceI = procPatch.start() + pfI;
                        label bfI = faceI - mesh_.nInternalFaces();

                        if
                        (
                            bfZoneID[bfI] == zone1ID
                         && coupleZoneID [bfI] == zone2ID
                        )
                        {
                            reportFaces_.set(faceI, 1);
                        }
                        else if
                        (
                            bfZoneID[bfI] == zone2ID
                         && coupleZoneID [bfI] == zone1ID
                        )
                        {
                            reportFaces_.set(faceI, 1);
                        }
                    }
                }
            }
        }
    }
    else
    {
        WarningInFunction
            << "cellZones " << zone1name << " or " << zone2name
            << " of cellZonePair could not be found"
            << " and will not be added to " << name()
            << " surface report." << endl;
    }
}


void functionObjects::surfaceReport::markPatchFaces(const dictionary& sD)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    labelHashSet patchSet;
    if (word(sD.lookup("type")) == "patch")
    {
        word patchName(sD.lookup("name"));
        label patchID = mesh_.boundaryMesh().findPatchID(patchName);
        if (patchID != -1)
        {
            patchSet.insert(patchID);
        }
        else
        {
            WarningInFunction
                << "Patch " << patchName << " could not be found"
                << " and will not be added to " << name()
                << " surface report." << endl;
        }
    }
    else
    {
        patchSet = pbm.patchSet(wordReList(sD.lookup("names")), false, true);
        if (patchSet.size() == 0)
        {
            WarningInFunction
                << "Patches " << sD.lookup("names") << " could not be found"
                << " and will not be added to " << name()
                << " surface report." << endl;
        }
    }

    forAllConstIter(labelHashSet, patchSet, iter)
    {
        label patchi = iter.key();
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        if (isA<indirectPolyPatch>(pp))
        {
            const indirectPolyPatch& idpp =
                refCast<const indirectPolyPatch>(pp);
            const labelList& addressing = idpp.addressing();
            forAll(addressing, ai)
            {
                label facei = addressing[ai];
                reportFaces_.set(facei, 1);
                if (facei < mesh_.nInternalFaces())
                {
                    reportContainsInternalFaces_ = true;
                }
            }
        }
        else
        {
            forAll(pp, pfI)
            {
                const label fI = pp.localToGlobal(pfI);
                reportFaces_.set(fI, 1);
                if (fI < mesh_.nInternalFaces())
                {
                    // Can happen for an indirect patch
                    reportContainsInternalFaces_ = true;
                }
            }
        }
    }
}


void functionObjects::surfaceReport::markPlaneFaces(const dictionary& sD)
{
    reportContainsInternalFaces_ = true;

    vector basePoint(sD.lookup("basePoint"));
    vector normal(sD.lookup("normal"));
    normal /= mag(normal);

    const labelList& neis = mesh_.neighbour();
    const labelList& owns = mesh_.faceOwner();

    scalarField toPlane((mesh_.C().primitiveField() - basePoint) & normal);

    forAll(toPlane, cI)
    {
        if (toPlane[cI] == 0)
        {
            toPlane[cI] = SMALL;
        }
    }

    forAll(neis, nI)
    {
        label neiCellI = neis[nI];
        label ownCellI = owns[nI];

        if (sign(toPlane[ownCellI]) != sign(toPlane[neiCellI]))
        {
            reportFaces_.set(nI, 1);
        }
    }

    // Processor boundary faces

    // Calculate coupled cell centres
    vectorField neiCc(mesh_.nFaces() - mesh_.nInternalFaces());

    for
    (
        label faceI = mesh_.nInternalFaces();
        faceI < mesh_.nFaces();
        faceI++
    )
    {
        neiCc[faceI - mesh_.nInternalFaces()] = mesh_.C()[owns[faceI]];
    }

    vectorField myCc = neiCc;
    syncTools::swapBoundaryFaceList(mesh_, neiCc);

    const polyPatchList& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            // Only add lower numbered proc face to avoid duplication
            if (procPatch.myProcNo() < procPatch.neighbProcNo())
            {
                forAll(procPatch, pfI)
                {
                    label faceI = procPatch.start() + pfI;
                    label bfI = faceI - mesh_.nInternalFaces();

                    scalar toPlaneThisProc = (myCc[bfI] - basePoint) & normal;
                    scalar toPlaneOtherProc = (neiCc[bfI] - basePoint) & normal;

                    if (toPlaneThisProc == 0)
                    {
                        toPlaneThisProc = SMALL;
                    }
                    if (toPlaneOtherProc == 0)
                    {
                        toPlaneOtherProc = SMALL;
                    }

                    if (sign(toPlaneThisProc) != sign(toPlaneOtherProc))
                    {
                        reportFaces_.set(faceI, 1);
                    }
                }
            }
        }
    }
}


void functionObjects::surfaceReport::setDirection()
{
    // Set the surfaceDirection to the normal of the first plane, otherwise
    // normal of the first face in the report

    surfaceDirectionPtr_.clear();

    PtrList<dictionary> surfDicts(dict_.lookup("surfaces"));
    forAll(surfDicts, dI)
    {
        const dictionary& sD = surfDicts[dI];

        if (word(sD.lookup("type")) == "plane")
        {
            vector basePoint(sD.lookup("basePoint"));
            vector normal(sD.lookup("normal"));
            normal /= mag(normal);
            surfaceDirectionPtr_.reset(new vector(normal));

            break;
        }
    }

    if (!surfaceDirectionPtr_.valid())
    {
        surfaceDirectionPtr_.reset(new vector(vector::zero));

        forAll(reportFaces_, rfI)
        {
            if (reportFaces_.get(rfI))
            {
                if (rfI < mesh_.nInternalFaces())
                {
                    surfaceDirectionPtr_.reset(new vector(mesh_.Sf()[rfI]));
                    break;
                }
                else
                {
                    label patchIndex = mesh_.boundaryMesh().whichPatch(rfI);
                    const polyPatch& pp = mesh_.boundaryMesh()[patchIndex];

                    label prfI = rfI - pp.start();

                    if (isA<nonConformalOrigPolyPatch>(pp))
                    {
                        const label polyBFacei = rfI - mesh_.nInternalFaces();
                        const labelUList& patches =
                            mesh_.polyBFacePatches()[polyBFacei];
                        const labelUList& patchFaces =
                            mesh_.polyBFacePatchFaces()[polyBFacei];

                        forAll(patches, patchi)
                        {
                            if (patches[patchi] != patchIndex)
                            {
                                patchIndex = patches[patchi];
                                prfI = patchFaces[patchi];

                                break;
                            }
                        }
                    }

                    surfaceDirectionPtr_.reset
                    (
                        new vector
                        (
                            mesh_.Sf().boundaryField()[patchIndex][prfI]
                        )
                    );

                    break;
                }
            }
        }

        // reduce(surfaceDirectionPtr_(), maxOp<vector>());
        reduceMaxVectorMag(surfaceDirectionPtr_());
        scalar magSurfaceDirectionPtr = mag(surfaceDirectionPtr_());
        if (magSurfaceDirectionPtr > SMALL)
        {
            surfaceDirectionPtr_() /= magSurfaceDirectionPtr;
        }
    }
}


void functionObjects::surfaceReport::updateFaces(const dictionary& dict)
{
    // Select faces
    markFaces(dict);

    // Unselect faces
    unmarkFaces(dict);

    // Calculate reporting direction
    setDirection();

    // Calculate surface weight
    calcWeighting();
}


void functionObjects::surfaceReport::markFaces(const dictionary& dict)
{
    PtrList<dictionary> surfDicts(dict.lookup("surfaces"));

    // Alternative interface for multiple patches
    wordReList patchSet
    (
        wordReList(dict.lookupOrDefault("patches", wordReList()))
    );

    forAll(patchSet, i)
    {
        dictionary entry;
        entry.add("type", "patch");
        entry.add("name", patchSet[i]);
        surfDicts.append(new dictionary(entry));
    }

    if (surfDicts.size() == 0)
    {
        WarningInFunction
            << "Deactivating " << name() << " as no surfaces entry present"
            << endl;
    }

    reportContainsInternalFaces_ = false;
    reportFaces_ = PackedList<2>(mesh_.nFaces(), 0);

    forAll(surfDicts, dI)
    {
        const dictionary& sD = surfDicts[dI];

        if (word(sD.lookup("type")) == "faceZone")
        {
            markFaceZone(sD);
        }
        else if (word(sD.lookup("type")) == "cellZonePair")
        {
            markCellZonePair(sD);
        }
        else if
        (
            word(sD.lookup("type")) == "patch"
         || word(sD.lookup("type")) == "patches"
        )
        {
            markPatchFaces(sD);
        }
        else if (word(sD.lookup("type")) == "plane")
        {
            markPlaneFaces(sD);
        }
        else
        {
            FatalErrorInFunction
                << "Invalid surface indicator for surface report: "
                << word(sD.lookup("type"))
                << nl << "Valid options are "
                << "'cellZonePair', 'faceZone', 'patch' and 'plane'"
                << exit(FatalError);
        }
    }

    reportContainsInternalFaces_ = returnReduce((bool) reportContainsInternalFaces_, orOp<bool>());
}


void functionObjects::surfaceReport::unmarkFaces(const dictionary& dict)
{
    // Unset non-master processor faces
    PackedBoolList isMasterFaces(syncTools::getMasterFaces(mesh_));
    forAll(mesh_.boundaryMesh(), pI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[pI];

        if (pp.coupled())
        {
            label pstart = pp.start();

            forAll(pp, pfI)
            {
                label faceI = pstart + pfI;

                if (reportFaces_.get(faceI) && !isMasterFaces[faceI])
                {
                    reportFaces_.set(faceI, 0);
                }
            }
        }
    }

    // Constrain face selection using bounding boxes
    PtrList<dictionary> conDicts(dict.lookup("constraints"));

    forAll(conDicts, dI)
    {
        const dictionary& cD = conDicts[dI];

        if (word(cD.lookup("type")) == "boundBox")
        {
            vector bbmin(cD.lookup("min"));
            vector bbmax(cD.lookup("max"));

            boundBox bb(bbmin, bbmax);

            Switch keepInside(cD.lookupOrDefault<Switch>("keepInside", true));

            if (keepInside)
            {
                forAll(mesh_.faceCentres(), fcI)
                {
                    if (!bb.contains(mesh_.faceCentres()[fcI]))
                    {
                        reportFaces_.set(fcI, 0);
                    }
                }
            }
            else
            {
                forAll(mesh_.faceCentres(), fcI)
                {
                    if (!bb.contains(mesh_.faceCentres()[fcI]))
                    {
                        reportFaces_.set(fcI, 0);
                    }
                }
            }
        }
        else if (word(cD.lookup("type")) == "sphere")
        {
            vector sc(cD.lookup("centre"));
            scalar sr(readScalar(cD.lookup("radius")));

            Switch keepInside(cD.lookupOrDefault<Switch>("keepInside", true));

            if (keepInside)
            {
                forAll(mesh_.faceCentres(), fcI)
                {
                    if (mag(mesh_.faceCentres()[fcI] - sc) > sr)
                    {
                        reportFaces_.set(fcI, 0);
                    }
                }
            }
            else
            {
                forAll(mesh_.faceCentres(), fcI)
                {
                    if (mag(mesh_.faceCentres()[fcI] - sc) <= sr)
                    {
                        reportFaces_.set(fcI, 0);
                    }
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Invalid constraint type for surface report: "
                << word(cD.lookup("type"))
                << nl << "Valid options are: "
                << "'boundBox'" << exit(FatalError);
        }
    }
}


void functionObjects::surfaceReport::calcWeighting()
{
    if (fluxWeighting_)
    {
        totalFlux();
    }

    if (!fluxWeighting_ || backFlowReport_)
    {
        totalSurface_= 0;

        const surfaceScalarField& magSf = mesh_.magSf();

        forAll(magSf, faceI)
        {
            if (reportFaces_.get(faceI))
            {
                totalSurface_ += magSf[faceI];
            }
        }

        forAll(magSf.boundaryField(), pI)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[pI];
            const label pstart = pp.start();

            if (!isA<nonConformalOrigPolyPatch>(pp))
            {
                forAll(magSf.boundaryField()[pI], pfI)
                {
                    if (reportFaces_.get(pstart + pfI))
                    {
                        totalSurface_ += magSf.boundaryField()[pI][pfI];
                    }
                }
            }
            else
            {
                // Add surface area from non-conformal patches
                forAll(magSf.boundaryField()[pI], pfI)
                {
                    if (reportFaces_.get(pstart + pfI))
                    {
                        const label polyBFacei =
                            pstart + pfI - mesh_.nInternalFaces();
                        const labelUList patches =
                            mesh_.polyBFacePatches()[polyBFacei];
                        const labelUList patchFaces =
                            mesh_.polyBFacePatchFaces()[polyBFacei];

                        forAll(patches, i)
                        {
                            totalSurface_ +=
                                magSf.boundaryField()
                                [patches[i]][patchFaces[i]];
                        }
                    }
                }
            }
        }

        reduce(totalSurface_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

functionObjects::surfaceReport::surfaceReport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),
    homogeneity_(false),
    fluxWeighting_(true),
    directionAware_(false),
    addDirectionMeanOutput_(false),
    fluxName_("phi"),
    fields_(0),
    surfaceDirectionPtr_(nullptr),
    reportFaces_(mesh_.nFaces(), 0),
    reportContainsInternalFaces_(false),
    compressible_(false),
    backFlowReport_(false),
    backFlowArea_(0.0),
    backFlowPer_(0.0),
    totalSurface_(GREAT),
    totalFlux_(GREAT),
    relativeFlux_(0.0),
    data_(),
    convergenceChecks_(),
    haveIndirectPatches_(false),
    meshChanged_(false),
    meshMoved_(false),
    fileWriter_
    (
        mesh_.time(),
        fileName
        (
            outputFileDir() + "/" + name + "/" +
            Time::timeName
            (
                mesh_.time().timeToUserTime(mesh_.time().startTime().value())
            )
        ),
        typeName,
        dict
    )
{
    // Check for indirect patches
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvP = mesh_.boundary()[patchi];
        if (isA<indirectPolyPatch>(fvP.patch()))
        {
            haveIndirectPatches_ = true;
            break;
        }
    }

    read(dict);
    writeFileHeader();

    // The legacy surfaceReport format requires that columns are
    // space-separated, but the header is tab-separated. This is all the
    // special handling that's required.
    if (fileWriter_.isDatFile())
    {
        fileWriter_.setDelimiter(token::SPACE);
        fileWriter_.setFixedWidth(false);
        fileWriter_.setNumberOfSpacesAfterDelimiter(0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

functionObjects::surfaceReport::~surfaceReport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool functionObjects::surfaceReport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    backFlowReport_ =
        dict.lookupOrDefault<Switch>("backFlowReport", "false");

    // Determine area/flux weighting.
    // Do this first because a lack of phi field deactivates the functionObject
    if (dict.found("weighting"))
    {
        word weighting = dict.lookup("weighting");

        if (weighting == "area")
        {
            fluxWeighting_ = false;
        }
        else if (weighting == "flux")
        {
            fluxWeighting_ = true;
            directionAware_ =
                dict.lookupOrDefault<bool>("directionAware", false);
            addDirectionMeanOutput_ =
                dict.lookupOrDefault<bool>("addDirectionMeanOutput", false);

            if (addDirectionMeanOutput_)
            {
                directionAware_ = false;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Invalid weighting factor for surface report: "
                << weighting
                << nl << "Valid weighting options for surfaceReport are "
                << "'area' and 'flux'" << exit(FatalError);
        }
    }

    const dictionary* subDictPtr = dict.subDictPtr("swirlDict");

    swirl_ = false;
    if (subDictPtr)
    {
        swirl_ = true;
        swirlUnitNorm_ = subDictPtr->lookup("swirlDir");
        swirlCentre_ = subDictPtr->lookup("swirlCentre");
    }

    const dictionary* tumbleDictPtr = dict.subDictPtr("tumbleDict");

    tumble_ = false;
    if (tumbleDictPtr)
    {
        tumble_ = true;
        tumbleDir_ = tumbleDictPtr->lookup("zDir");
        tumbleOrig_ = tumbleDictPtr->lookup("origin");
        rhoName_ = tumbleDictPtr->lookupOrDefault<word>("rhoName", "rhoInf");
        rhoRef_ = tumbleDictPtr->lookupOrDefault<scalar>("rhoRef", 1.205);
        Ls_ = tumbleDictPtr->lookupOrDefault<scalar>("Ls", 0.0932);
        massAveraged_ =
            tumbleDictPtr->lookupOrDefault<bool>("massAveraged", false);
    }

    // Determine compressible/incompressible
    {
        fluxName_ = dict.lookupOrDefault<word>("fluxName", "phi");

        if
        (
            !foundObject<surfaceScalarField>(fluxName_)
        )
        {
            WarningInFunction
                << fluxName_ << " could not be found in the database, "
                << "deactivating."
                << endl;
        }
        else if
        (
            lookupObject<surfaceScalarField>(fluxName_).dimensions()
         == dimDensity*dimVelocity*dimArea
        )
        {
            compressible_ = true;
        }
    }

    homogeneity_= false;

    if (dict.found("homogeneity"))
    {
        homogeneity_ = readBool(dict.lookup("homogeneity"));
    }

    // Field names
    {
        wordList tfns = dict.lookup("fields");
        fields_.setSize(tfns.size());
        fields_ = tfns;
        data_.setSize(tfns.size());

        forAll(data_, dI)
        {
            data_[dI] = -1;
        }
    }

    const dictionary* convSubDictPtr = dict.subDictPtr("convergence");

    if (convSubDictPtr)
    {
        label nConvChecks = 0;
        forAll(fields_, fieldI)
        {
            word fieldName = fields_[fieldI];
            if (convSubDictPtr->found(fieldName))
            {
                nConvChecks++;
            }
        }
        convergenceChecks_.setSize(nConvChecks);

        nConvChecks = 0;
        forAll(fields_, fieldI)
        {
            word fieldName = fields_[fieldI];
            if (convSubDictPtr->found(fieldName))
            {
                const dictionary convDict(convSubDictPtr->subDict(fieldName));

                convergenceChecks_.set
                (
                    nConvChecks,
                    new convergenceTermination(obr_, convDict)
                );
                nConvChecks++;
            }
        }
    }

    // Select faces
    updateFaces(dict);

    return true;
}


bool Foam::functionObjects::surfaceReport::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    // Main driver function
    calculate();

    // Write to screen

    // Net flux
    if (compressible_)
    {
        Log << "massflux [kg/s] = ";
    }
    else
    {
        Log << "volumeflux [m^3/s] = ";
    }
    Log << relativeFlux_ << endl;

    if (backFlowReport_)
    {
        if (totalSurface_ > 0)
        {
            backFlowPer_ = backFlowArea_/totalSurface_*100;
            Log << "backflow in " << backFlowArea_ << "[m^2] out of "
                << totalSurface_ << " [m^2] - " << backFlowPer_ << "%."
                << endl;
        }
        else
        {
            backFlowPer_ = 0.0;
            Log << "Total surface is zero. BackFlow not applicable." << endl;
        }
    }

    forAll(fields_, fI)
    {
        if (homogeneity_)
        {
            if (addDirectionMeanOutput_)
            {
                Log <<  fields_[fI]
                    << ": min " << data_[fI][0]
                    << ", max " << data_[fI][1]
                    << ", mean " << data_[fI][2]
                    << ", directionMean " << data_[fI][4]
                    << ", homogeneity " << data_[fI][3] << endl;
            }
            else
            {
                Log <<  fields_[fI]
                    << ": min " << data_[fI][0]
                    << ", max " << data_[fI][1]
                    << ", mean " << data_[fI][2]
                    << ", homogeneity " << data_[fI][3] << endl;
            }
        }
        else
        {
            if (addDirectionMeanOutput_)
            {
                Log <<  fields_[fI]
                    << ": min " << data_[fI][0]
                    << ", max " << data_[fI][1]
                    << ", mean " << data_[fI][2]
                    << ", directionMean " << data_[fI][4]
                    << ", stdDev " << data_[fI][3] << endl;
            }
            else
            {
                Log <<  fields_[fI]
                    << ": min " << data_[fI][0]
                    << ", max " << data_[fI][1]
                    << ", mean " << data_[fI][2]
                    << ", stdDev " << data_[fI][3] << endl;
            }
        }
    }

    if (swirl_)
    {
        Log <<  "Calculated Swirl " << swirlData_ << endl;
    }

    if (tumble_)
    {
        Log <<  "Calculated Tumble " << tumbleValue_ << endl;
    }

    fileWriter_.writeTime();
    fileWriter_.writeDelimited(relativeFlux_);

    forAll(fields_, fI)
    {
        if (addDirectionMeanOutput_)
        {
            fileWriter_.writeDelimited
            (
                data_[fI][0],
                data_[fI][1],
                data_[fI][2],
                data_[fI][4],
                data_[fI][3]
            );
        }
        else
        {
            fileWriter_.writeDelimited
            (
                data_[fI][0],
                data_[fI][1],
                data_[fI][2],
                data_[fI][3]
            );
        }

        word nameStr('(' + fields_[fI] + ')');
        this->setResult("min" + nameStr, data_[fI][0]);
        this->setResult("max" + nameStr, data_[fI][1]);
        this->setResult("mean" + nameStr, data_[fI][2]);

        if (addDirectionMeanOutput_)
        {
            this->setResult("directionMean" + nameStr, data_[fI][4]);
        }

        if (homogeneity_)
        {
            this->setResult("homogeneity" + nameStr, data_[fI][3]);
        }
        else
        {
            this->setResult("stdDev" + nameStr, data_[fI][3]);
        }
    }

    if (swirl_)
    {
        fileWriter_.writeDelimited(swirlData_);
    }
    if (tumble_)
    {
        fileWriter_.writeDelimited(tumbleValue_);
    }
    if (backFlowReport_)
    {
        fileWriter_.writeDelimited(backFlowPer_);
    }

    fileWriter_.endLine();
    Log << endl;

    return true;
}


void functionObjects::surfaceReport::calculate()
{
    if (haveIndirectPatches_ || meshChanged_ || meshMoved_)
    {
        // Topo changes or mesh movement: re-mark the surface
        updateFaces(dict_);
        meshChanged_ = false;
        meshMoved_ = false;
    }
    else if (fluxWeighting_)
    {
        totalFlux();
    }

    relativeFlux();

    if (backFlowReport_)
    {
        backFlowArea();
    }

    if (swirl_)
    {
        if (foundObject<volVectorField>("U"))
        {
            calcSwirlData(lookupObject<volVectorField>("U"));
        }
        else
        {
            WarningInFunction
                << "Field U could not be found in the object registry."
                << " No swirl report compiled." << endl;
        }
    }

    if (tumble_)
    {
        if (foundObject<volVectorField>("U"))
        {
            calcTumble(lookupObject<volVectorField>("U"));
        }
        else
        {
            WarningInFunction
                << "Field U could not be found in the object registry."
                << " No tumble report compiled." << endl;
        }
    }

    forAll(fields_, fI)
    {
        if (foundObject<volScalarField>(fields_[fI]))
        {
            calcScalarData
            (
                lookupObject<volScalarField>(fields_[fI]),
                data_[fI]
            );
        }
        else if (foundObject<volVectorField>(fields_[fI]))
        {
            calcNonScalarData
            (
                lookupObject<volVectorField>(fields_[fI]),
                data_[fI]
            );
        }
        else if (foundObject<volTensorField>(fields_[fI]))
        {
            calcNonScalarData
            (
                lookupObject<volTensorField>(fields_[fI]),
                data_[fI]
            );
        }
        else if (foundObject<volSymmTensorField>(fields_[fI]))
        {
            calcNonScalarData
            (
                lookupObject<volSymmTensorField>(fields_[fI]),
                data_[fI]
            );
        }
        else if (foundObject<volSphericalTensorField>(fields_[fI]))
        {
            calcNonScalarData
            (
                lookupObject<volSphericalTensorField>(fields_[fI]),
                data_[fI]
            );
        }
        else
        {
            WarningInFunction
                << "Field " << fields_[fI]
                << " could not be found in the object registry."
                << " No surface report compiled." << endl;
        }
    }

    forAll(fields_, fI)
    {
       word fieldName = fields_[fI];
       forAll(convergenceChecks_, checkI)
       {
           convergenceTermination& check = convergenceChecks_[checkI];

           if (check.fieldName() == fieldName)
           {
               if (check.calculate(data_[fI][2]))
               {
                   Log << "Convergence obatained for field : " << fieldName
                       << " in function object : " << name() << endl;
               }
           }
       }
   }
}


bool functionObjects::surfaceReport::write()
{
    return true;
}


void functionObjects::surfaceReport::updateMesh(const mapPolyMesh& mpm)
{
    fvMeshFunctionObject::updateMesh(mpm);
    meshChanged_ = true;
}


void functionObjects::surfaceReport::movePoints(const polyMesh& mesh)
{
    fvMeshFunctionObject::movePoints(mesh);
    meshMoved_ = true;
}


// ************************************************************************* //
