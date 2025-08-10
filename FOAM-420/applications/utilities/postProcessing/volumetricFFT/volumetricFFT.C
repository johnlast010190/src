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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2016 Esi Ltd.


 Application
 volumetricFFT.C

 Description
 Fourier transform analysis of computational domain

 \*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "switchParallel.H"
#include "fvMesh/fvMesh.H"
#include "noiseParameters.H"
#include "primitives/strings/wordRes/wordReListMatcher.H"

using namespace Foam ;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fftw3.h"

void readPressure
(
    const fvMesh& mesh,
    const noiseParameters& noiseParam,
    const instantList& readTimes,
    const label blockI,
    const scalar rho,
    const labelHashSet& calculatedPatches,
    scalarField& allTimeInternalData,
    scalarField& allTimePatchData,
    Time& runTime
)
{
    label timestep = -1 ;

    word pName = noiseParam.pName();
    bool internal = noiseParam.internal();
    label blockSize = noiseParam.blockSize();
    bool hanning = noiseParam.hanning();
    scalar hanningCoeff = noiseParam.hanningCoeff();
    label ntimesteps = noiseParam.numberOfTimeSteps();

    forAll(readTimes, timeI)
    {
        if ((timeI % 1000) == 0)
        {
            Pout<<"Read in pressure data for time: "<< timeI <<endl;
        }

        runTime.setTime(readTimes[timeI],timeI);
        ++timestep;

        if (blockSize > 0)
        {
            volScalarField p
            (
                IOobject
                (
                    pName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            if (blockI == 0)
            {
                label allTimeIndex = 0 ;
                label index = 0 ;

                forAll(p.boundaryField( ),patchI)
                {
                    if (!calculatedPatches.found(patchI))
                    {
                        continue;
                    }

                    for
                    (
                        fvPatchScalarField::iterator it =
                            p.boundaryFieldRef()[patchI].begin();
                        it != p.boundaryField()[patchI].end();
                        it++
                    )
                    {
                        allTimeIndex =
                            2*(timestep+index*ntimesteps);
                        if (hanning)
                        {
                            scalar theta =
                                Foam::constant::mathematical::twoPi
                                *timestep/(ntimesteps-1);

                            allTimePatchData[allTimeIndex] =
                                (*it)*rho*hanningCoeff*0.5
                                *(1-Foam::cos(theta));
                            allTimePatchData[allTimeIndex+1] = 0.0;
                        }
                        else
                        {
                            allTimePatchData[allTimeIndex] =
                                (*it)*rho;
                            allTimePatchData[allTimeIndex+1] = 0.0;
                        }
                        index++;
                    }
                }
            }

            if (internal)
            {
                label index = 0;
                label allTimeIndex = 0;

                for
                (
                    volScalarField::iterator it = p.begin()+(blockI*blockSize);
                    it!= p.end() && ((index%blockSize) || index == 0);
                    it++
                )
                {
                    allTimeIndex = 2 * ( timestep + index * ntimesteps );
                    if (hanning)
                    {
                        scalar theta = Foam::constant::mathematical::twoPi
                            *timestep/(ntimesteps-1);
                        allTimeInternalData[allTimeIndex] = (*it)*rho
                            *hanningCoeff*0.5*(1-Foam::cos(theta));
                        allTimeInternalData[allTimeIndex+1] = 0.0;
                    }
                    else
                    {
                        allTimeInternalData[allTimeIndex] = (*it)*rho;
                        allTimeInternalData[allTimeIndex+1] =0.0;
                    }
                    index++ ;
                }
            }
        }
    }
}


void calculateDFT
(
    const noiseParameters& noiseParam,
    const label blockI,
    const label nBoundaryFaces,
    fftw_complex* fx,
    fftw_complex* fxSurf
)
{
    label ntimesteps = noiseParam.numberOfTimeSteps();
    bool internal = noiseParam.internal();
    label blockSize = noiseParam.blockSize();

    if (blockSize > 0)
    {
        fftw_plan planSurf;
        fftw_plan plan;

        if (blockI == 0)
        {
            for (int faceIndex = 0 ; faceIndex < nBoundaryFaces ; faceIndex++)
            {
                planSurf = fftw_plan_dft_1d
                (
                    ntimesteps,
                    fxSurf + faceIndex * ntimesteps,
                    fxSurf + faceIndex * ntimesteps,
                    FFTW_FORWARD ,
                    FFTW_ESIIMATE
                 );
                fftw_execute(planSurf);
                fftw_destroy_plan(planSurf);
            }
        }

        if (internal)
        {
            for
            (
                int cellIndex = 0;
                cellIndex == 0 || (cellIndex%blockSize);
                cellIndex++
            )
            {
                plan = fftw_plan_dft_1d
                (
                    ntimesteps,
                    fx + cellIndex * ntimesteps,
                    fx + cellIndex * ntimesteps,
                    FFTW_FORWARD,
                    FFTW_ESIIMATE
                  );
                fftw_execute(plan);
                fftw_destroy_plan(plan);
            }
        }
    }
}


void calculateOneThirds
(
    const fvMesh& mesh,
    const noiseParameters& noiseParam,
    const Time& runTime,
    const fftw_complex* fx,
    const fftw_complex* fxSurf,
    const label blockI,
    const label ncells,
    const labelHashSet& calculatedPatches,
    const instantList& readTimes,
    const scalar deltaF,
    DynamicList<scalar>& octfm,
    DynamicList<DynamicList<scalar>>& freqPerOneThirdOctave
)
{
    scalar P0 = 2.0e-5;

    label ntimesteps = noiseParam.numberOfTimeSteps();
    bool internal = noiseParam.internal();
    label blockSize = noiseParam.blockSize();

    label timeI = 0 ;
    scalar realcoeff = ntimesteps/2;
    scalar imgcoeff = -realcoeff;

    if (blockI == 0)
    {
        forAll(freqPerOneThirdOctave, fm)
        {
            volScalarField fft
            (
                IOobject
                (
                    "octave_"+name(octfm(fm)),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "octave_"+name(octfm(fm)),
                    dimless,
                    0.0
                 )
            );
            fft.write();
        }
    }

    forAll(freqPerOneThirdOctave, fm)
    {
        volScalarField fft
        (
            IOobject
            (
                "octave_"+name(octfm(fm)),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
             ),
             mesh
          );

        forAll(freqPerOneThirdOctave(fm), f)
        {
            timeI = freqPerOneThirdOctave(fm)(f);

            if (timeI == 0 || timeI == ntimesteps)
            {
                realcoeff = ntimesteps;
            }
            else
            {
                realcoeff = ntimesteps/2;
            }

            label faceIndex = -1;

            if (blockSize > 0)
            {
                if (blockI == 0)
                {
                    volScalarField::Boundary& fftbf = fft.boundaryFieldRef();
                    forAll(fftbf, patchI)
                    {
                        if (!calculatedPatches.found(patchI))
                        {
                            continue;
                        }
                        forAll(fftbf[patchI], faceI)
                        {
                            faceIndex++ ;

                            scalar a1 = Foam::sqr
                            (
                                fxSurf[timeI+faceIndex*ntimesteps][0]
                                /realcoeff
                            );

                            scalar a2 = Foam::sqr
                            (
                                fxSurf[timeI+faceIndex*ntimesteps][1]
                                /imgcoeff
                            );

                            if (Foam::sqrt(a1+a2) > 0)
                            {
                                fftbf[patchI][faceI] +=
                                    Foam::pow
                                    (
                                        10.0,
                                        2.0*Foam::log10
                                        (Foam::sqrt(a1+a2)/P0)
                                     );
                            }
                        }
                    }
                }

                if (internal)
                {
                    for
                    (
                        int cellIndex = 0;
                        (
                            cellIndex == 0
                            || (cellIndex%blockSize)
                        ) && cellIndex+blockI*blockSize< ncells;
                        ++cellIndex
                     )
                    {
                        scalar a1 = Foam::sqr
                        (
                            fx[timeI+cellIndex*ntimesteps][0]
                            /realcoeff
                        );
                        scalar a2 = Foam::sqr
                        (
                            fx[timeI+cellIndex*ntimesteps][1]
                            /imgcoeff
                         );

                        if (Foam::sqrt(a1+a2) > 0)
                        {
                            fft[cellIndex + blockI * blockSize] =
                                fft[cellIndex+blockI*blockSize]
                                + Foam::pow
                                (
                                    10.0 ,
                                    2.0*Foam::log10(Foam::sqrt(a1+a2)/P0)
                                );
                        }
                    }
                }
            }
        }

        if (blockSize > 0)
        {
            if (blockI == 0)
            {
                volScalarField::Boundary& fftbf = fft.boundaryFieldRef();
                forAll(fftbf, patchI)
                {
                    if (!calculatedPatches.found(patchI))
                    {
                        continue;
                    }
                    forAll(fftbf[patchI], faceI)
                    {
                        if (fftbf[patchI][faceI] > 0)
                        {
                            fftbf[patchI][faceI] =
                                10.0*Foam::log10(fftbf[patchI][faceI]);
                        }
                    }
                }
            }
            if (internal)
            {
                for
                (
                    int cellIndex = 0;
                    (cellIndex == 0 || ( cellIndex % blockSize ))
                        && cellIndex+blockI*blockSize < ncells;
                    ++cellIndex
                )
                {
                    if (fft[cellIndex+blockI*blockSize]>0)
                    {
                        fft[cellIndex+blockI*blockSize] = 10.0*Foam::log10
                            (fft[cellIndex+blockI*blockSize]);
                    }
                }
            }
        }
        fft.write();
    }
}


void calculateAllFrequencies
(
    const fvMesh& mesh,
    const noiseParameters& noiseParam,
    const instantList& readTimes,
    const fftw_complex* fx,
    const fftw_complex* fxSurf,
    const label blockI,
    const label ncells,
    const labelList& patchIDs,
    Time& runTime
)
{
    scalar P0 = 2.0e-5;
    label timestep = -1 ;

    label ntimesteps = noiseParam.numberOfTimeSteps();
    bool internal = noiseParam.internal();
    label blockSize = noiseParam.blockSize();

    labelHashSet calculatedPatches(patchIDs);

    forAll(readTimes , timeI)
    {
        runTime.setTime(readTimes[timeI], timeI);
        ++timestep;

        scalar realcoeff = ntimesteps/2;
        scalar imgcoeff = -realcoeff;

        autoPtr<volScalarField> fftdBPtr;

        if (patchIDs.size() == 0)
        {
            fftdBPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "fftdB",
                        runTime.timeName( ),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                      ),
                    mesh,
                    dimensionedScalar(0.0),
                    zeroGradientFvPatchScalarField::typeName
                 )
             );
        }
        else
        {
            fftdBPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "fftdB",
                        runTime.timeName( ),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(0.0)
                 )
             );
        }
        volScalarField& fftdB = fftdBPtr();

        label faceIndex = -1 ;

        if (blockSize > 0)
        {
            if (blockI == 0)
            {
                volScalarField::Boundary& fftdBbf = fftdB.boundaryFieldRef();
                forAll(fftdBbf, patchI)
                {
                    if (!calculatedPatches.found(patchI))
                    {
                        continue;
                    }
                    forAll(fftdBbf[patchI], faceI)
                    {
                        faceIndex++;

                        scalar a1 = Foam::sqr
                        (
                            fxSurf
                            [timestep+faceIndex*ntimesteps][0]
                            /realcoeff
                         );
                        scalar a2 = Foam::sqr
                        (
                            fxSurf
                            [timestep+faceIndex*ntimesteps][1]
                            /imgcoeff
                         );

                        if (Foam::sqrt(a1+a2) > 0)
                        {
                            fftdBbf[patchI][faceI]
                                = 20.0*Foam::log10(Foam::sqrt(a1+a2)/P0);
                        }
                    }
                }
            }

            if (internal)
            {
                for
                (
                    int cellIndex = 0;
                    (cellIndex == 0 || (cellIndex%blockSize))
                        && cellIndex+blockI*blockSize < ncells;
                    ++cellIndex
                )
                {
                    if
                    (
                        timestep == 0 || timestep == ntimesteps
                    )
                    {
                        realcoeff = ntimesteps ;
                    }
                    else
                    {
                        realcoeff = ntimesteps / 2 ;
                    }

                    label index = timestep+cellIndex* ntimesteps;
                    scalar a1 =
                    (
                        (fx[index][0]/realcoeff)*(fx[index][0]/realcoeff)
                        +(fx[index][1]/imgcoeff)*(fx[index][1]/imgcoeff)
                    );

                    if (Foam::sqrt(a1) > 0)
                    {
                        fftdB[cellIndex + blockI * blockSize] = 20.0
                            * Foam::log10(Foam::sqrt(a1)/P0);
                    }
                }
            }
        }
        fftdB.write( ) ;
    }
}


void calculateReverseFFT
(
    const fvMesh& mesh,
    const noiseParameters& noiseParam,
    const instantList& readTimes,
    const label blockI,
    const scalar rho,
    const label nBoundaryFaces,
    const label ncells,
    const labelList& patchIDs,
    const scalar deltaF,
    fftw_complex* fx,
    fftw_complex* fxSurf,
    Time& runTime
)
{
    fftw_plan planSurf;
    fftw_plan plan;

    label ntimesteps = noiseParam.numberOfTimeSteps();
    bool internal = noiseParam.internal();
    label blockSize = noiseParam.blockSize();
    bool hanning = noiseParam.hanning();
    scalar hanningCoeff = noiseParam.hanningCoeff();
    label reverseFFTRatio = noiseParam.reverseFFTRatio();
    const vector2DField& filterBands = noiseParam.filterBands();

    labelHashSet calculatedPatches(patchIDs);

    forAll(filterBands, filterI)
    {
        scalarField outSurfaceData(2*ntimesteps*nBoundaryFaces, 0.0);
        fftw_complex * outSurf = reinterpret_cast <fftw_complex*>
            (outSurfaceData.begin());

        scalarField outData(0);
        if (internal)
        {
            outData.setSize(2*blockSize*ntimesteps, 0.0);
        }
        fftw_complex * out = reinterpret_cast <fftw_complex*>
            (outData.begin());

        label count = 0;

        List<scalar> removedSurf(0);
        List<scalar> removed(0);

        if (blockSize > 0)
        {
            if (blockI == 0)
            {
                for
                (
                    int faceIndex = 0 ;
                    faceIndex < nBoundaryFaces ;
                    faceIndex++
                )
                {

                    for
                    (
                        int fnum = 0;
                        fnum < ntimesteps / 2;
                        fnum++
                    )
                    {
                        bool keep = false ;

                        if
                        (
                            fnum > int(filterBands[filterI].x()/deltaF)
                            && fnum < int(filterBands[filterI].y()/deltaF)
                        )
                        {
                            keep = true ;
                        }

                        if (!keep)
                        {
                            count += 4;
                        }
                    }
                }
            }
            removedSurf.setSize(count);

            count = 0;
            if (internal)
            {
                for
                (
                    int cellIndex = 0;
                    cellIndex == 0 || ( cellIndex % blockSize );
                    ++cellIndex
                )
                {
                    for
                    (
                        int fnum = 0;
                        fnum < ntimesteps / 2;
                        fnum++
                    )
                    {
                        bool keep = false ;

                        if
                        (
                            fnum > int(filterBands[filterI].x()/deltaF)
                            && fnum < int(filterBands[filterI].y()/deltaF)
                        )
                        {
                            keep = true ;
                        }

                        if (!keep)
                        {
                            count += 4;
                        }
                    }
                }
            }

            removed.setSize(count);

            count = 0;
            if (blockI == 0)
            {
                for
                (
                    int faceIndex = 0 ;
                    faceIndex < nBoundaryFaces ;
                    faceIndex++
                )
                {

                    for
                    (
                        int fnum = 0;
                        fnum < ntimesteps / 2;
                        fnum++
                    )
                    {
                        bool keep = false ;

                        if
                        (
                            fnum > int(filterBands[filterI].x()/deltaF)
                            && fnum < int(filterBands[filterI].y()/deltaF)
                        )
                        {
                            keep = true ;
                        }

                        if (!keep)
                        {
                            removedSurf[count++] =
                                fxSurf[faceIndex*ntimesteps+fnum][0];
                            fxSurf[faceIndex*ntimesteps+fnum][0] = 0.0;

                            removedSurf[count++] =
                                fxSurf[faceIndex*ntimesteps+fnum][1];
                            fxSurf[faceIndex*ntimesteps+fnum][1] = 0.0;

                            removedSurf[count++] =
                                fxSurf[(faceIndex+1)*ntimesteps-1-fnum][0];
                            fxSurf[(faceIndex+1)*ntimesteps-1-fnum][0]=0.0;

                            removedSurf[count++] =
                                fxSurf[(faceIndex+1)*ntimesteps-1-fnum][1];
                            fxSurf[(faceIndex+1)*ntimesteps-1-fnum][1]=0.0;
                        }
                    }

                    planSurf = fftw_plan_dft_1d
                    (
                        ntimesteps ,
                        fxSurf+faceIndex*ntimesteps,
                        outSurf+faceIndex*ntimesteps,
                        FFTW_BACKWARD,
                        FFTW_ESIIMATE
                    );
                    fftw_execute(planSurf);
                    fftw_destroy_plan(planSurf);
                }
            }

            count = 0;
            if (internal)
            {
                for
                (
                    int cellIndex = 0;
                    cellIndex == 0 || ( cellIndex % blockSize );
                    ++cellIndex
                )
                {
                    for
                    (
                        int fnum = 0;
                        fnum < ntimesteps / 2;
                        fnum++
                    )
                    {
                        bool keep = false ;

                        if
                        (
                            fnum > int(filterBands[filterI].x()/deltaF)
                            && fnum < int(filterBands[filterI].y()/deltaF)
                        )
                        {
                            keep = true ;
                        }

                        if (!keep)
                        {
                            removed[count++] =
                                fx[cellIndex*ntimesteps+fnum][0];
                            fx[cellIndex*ntimesteps+fnum][0]=0.0 ;

                            removed[count++] =
                                fx[cellIndex*ntimesteps+fnum][1];
                            fx[cellIndex*ntimesteps+fnum][1]=0.0 ;

                            removed[count++] =
                                fx[(cellIndex+1)*ntimesteps-1-fnum][0];
                            fx[(cellIndex+1)*ntimesteps-1-fnum][0] = 0.0;

                            removed[count++] =
                                fx[(cellIndex+1)*ntimesteps-1-fnum][1];
                            fx[(cellIndex+1)*ntimesteps-1-fnum][1] = 0.0 ;
                        }
                    }

                    plan = fftw_plan_dft_1d
                    (
                        ntimesteps,
                        fx+cellIndex*ntimesteps,
                        out+cellIndex*ntimesteps,
                        FFTW_BACKWARD,
                        FFTW_ESIIMATE
                    );
                    fftw_execute(plan);
                    fftw_destroy_plan(plan);
                }
            }
        }

        label timestep = -1 ;

        forAll(readTimes, timeI)
        {
            if (!(timeI % reverseFFTRatio))
            {
                runTime.setTime(readTimes[timeI],timeI);
                ++timestep;

                label filterBand = label
                (
                    (filterBands[filterI].x()
                     +filterBands[filterI].y())/2.
                );

                autoPtr<volScalarField> inversePtr;

                if (patchIDs.size() == 0)
                {
                    inversePtr.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "inverse_" + name(filterBand),
                                runTime.timeName( ) ,
                                mesh ,
                                IOobject::READ_IF_PRESENT ,
                                IOobject::AUTO_WRITE
                            ),
                            mesh,
                            dimensionedScalar
                            (
                                "inverse" ,
                                dimensionSet(1,-1,-2,0,0),
                                0.0
                            ),
                            zeroGradientFvPatchScalarField::typeName
                         )
                     );
                }
                else
                {
                    inversePtr.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "inverse_" + name(filterBand),
                                runTime.timeName( ) ,
                                mesh ,
                                IOobject::READ_IF_PRESENT ,
                                IOobject::AUTO_WRITE
                            ),
                            mesh,
                            dimensionedScalar
                            (
                                "inverse" ,
                                dimensionSet(1,-1,-2,0,0),
                                0.0
                             )
                         )
                     );
                }
                volScalarField& inverse = inversePtr();

                if (blockSize > 0)
                {
                    scalar hanningCorrection = 1.0/rho;
                    if (hanning)
                    {
                        scalar theta =
                            Foam::constant::mathematical::twoPi
                            *timestep/(ntimesteps - 1);

                        hanningCorrection =
                        (
                            rho*hanningCoeff*0.5*
                            (
                                1 - Foam::cos(theta)
                            )
                        ) ;
                    }

                    label faceIndex = -1 ;
                    if (blockI == 0)
                    {
                        volScalarField::Boundary& inversebf =
                            inverse.boundaryFieldRef();
                        forAll(inversebf, patchI)
                        {
                            if (!calculatedPatches.found(patchI))
                            {
                                continue;
                            }

                            forAll(inversebf[patchI], faceI)
                            {
                                faceIndex++ ;

                                if (hanningCorrection > 0)
                                {
                                    inversebf[patchI]
                                        [faceI] = outSurf
                                        [timestep+faceIndex*ntimesteps][0]/
                                        (
                                            ntimesteps
                                            *hanningCorrection
                                        );
                                }
                                else
                                {
                                    inversebf[patchI][faceI] = 0.0;
                                }
                            }
                        }
                    }

                    if (internal)
                    {
                        for
                        (
                            int cellIndex = 0;
                            (cellIndex == 0||(cellIndex%blockSize))
                                && cellIndex + blockI * blockSize < ncells;
                            ++cellIndex
                        )
                        {
                            if (hanningCorrection > 0)
                            {
                                inverse[cellIndex + blockI * blockSize] =
                                (
                                    out[timestep+cellIndex*ntimesteps][0]
                                    / ntimesteps
                                  )/ hanningCorrection ;
                            }
                            else
                            {
                                inverse[cellIndex+blockI*blockSize] =
                                    out[timestep+cellIndex*ntimesteps][0]
                                    / ntimesteps;
                            }
                        }
                    }
                }
                inverse.write( ) ;
            }
        }

        if (blockSize > 0)
        {
            //Reset filtered values
            count = 0;
            if (blockI == 0)
            {
                for
                (
                    int faceIndex = 0 ;
                    faceIndex < nBoundaryFaces ;
                    faceIndex++
                )
                {

                    for
                    (
                        int fnum = 0;
                        fnum < ntimesteps / 2;
                        fnum++
                    )
                    {
                        bool keep = false ;

                        if
                        (
                            fnum > int(filterBands[filterI].x()/deltaF)
                            && fnum < int(filterBands[filterI].y()/deltaF)
                        )
                        {
                            keep = true ;
                        }

                        if (!keep)
                        {
                            fxSurf[faceIndex*ntimesteps+fnum][0] =
                                removedSurf[count++];
                            fxSurf[faceIndex*ntimesteps+fnum][1] =
                                removedSurf[count++];
                            fxSurf[(faceIndex+1)*ntimesteps-1-fnum][0] =
                                removedSurf[count++];
                            fxSurf[(faceIndex+1)*ntimesteps-1-fnum][1] =
                                removedSurf[count++];
                        }
                    }
                }
            }

            count = 0;
            if (internal)
            {
                for
                (
                    int cellIndex = 0;
                    cellIndex == 0 || ( cellIndex % blockSize );
                    ++cellIndex
                )
                {
                    for
                    (
                        int fnum = 0;
                        fnum < ntimesteps / 2;
                        fnum++
                    )
                    {
                        bool keep = false ;

                        if
                        (
                            fnum > int(filterBands[filterI].x()/deltaF)
                            && fnum < int(filterBands[filterI].y()/deltaF)
                        )
                        {
                            keep = true ;
                        }

                        if (!keep)
                        {
                            fx[cellIndex*ntimesteps+fnum][0] =
                                removed[count++];
                            fx[cellIndex*ntimesteps+fnum][1] =
                                removed[count++];
                            fx[(cellIndex+1)*ntimesteps-1-fnum][0] =
                                removed[count++];
                            fx[(cellIndex+1)*ntimesteps-1-fnum][1] =
                                removed[count++];
                        }
                    }
                }
            }
        }
    }
}


void calculateOctaveFrequencies
(
    const noiseParameters& noiseParam,
    const instantList& readTimes,
    const scalar deltaF,
    DynamicList<scalar>& octfm,
    DynamicList<DynamicList<scalar>>& freqPerOneThirdOctave
)
{
    label ntimesteps = noiseParam.numberOfTimeSteps();
    scalar fLower = noiseParam.fLower();
    scalar fUpper = noiseParam.fUpper();

    DynamicList<scalar> octfl;
    DynamicList<scalar> octfu;
    scalarList allFrequencies(ntimesteps);

    scalar cubert2 = Foam::cbrt(2.0);
    scalar sixthrt2 = Foam::sqrt(cubert2);

    label pp = 0;

    while (pp == 0 || octfm(pp - 1) < fUpper)
    {
        octfm.append(fLower*pow(cubert2,pp));
        octfl.append(1.0/sixthrt2*octfm(pp));
        octfu.append(octfl(pp)*cubert2);
        pp++ ;
    }

    Info<<"Octave bands: index, mid, lower, higher"<<endl;
    forAll(octfm, i)
    {
        Info<<i+1<<" "<<octfm(i)<<" "<<octfl(i)
            <<" "<<octfu(i)<<endl;
    }

    freqPerOneThirdOctave.setSize(octfm.size());

    label octCount = 0 ;

    forAll(readTimes, timeI)
    {
        allFrequencies[timeI] = timeI*deltaF;
        if
        (
            allFrequencies[timeI] > octfl(0)
            && allFrequencies[timeI]
            < octfu(octfu.size()-1)
          )
        {
            if
            (
                allFrequencies[timeI] < octfu(octCount)
                && allFrequencies[timeI] >= octfl(octCount)
              )
            {
                freqPerOneThirdOctave(octCount).append(timeI);
            }
            else
            {
                octCount++ ;
                freqPerOneThirdOctave(octCount).append(timeI);
            }
        }
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "include/setRootCase.H"
    #include "include/createTime.H"

    instantList allTimes = timeSelector::select0( runTime , args ) ;

    #include "include/createMesh.H"

    IOdictionary noiseDict
    (
        IOobject
        (
            "volumetricFFTDict" ,
            runTime.system( ) ,
            mesh ,
            IOobject::MUST_READ ,
            IOobject::NO_WRITE
         )
    );

    noiseParameters noiseParam(noiseDict);

    if (noiseParam.binary())
    {
        runTime.setIOStream
        (
            IOstream::BINARY,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED
         );
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    DynamicList<label> patchIDs(patches.size());

    if (noiseParam.patches().size())
    {
        List<wordRe> patchRe(noiseParam.patches());
        wordReListMatcher patchesToWrite(patchRe);

        forAll(patches, i)
        {
            if (patchesToWrite.match(patches[i].name()))
            {
                patchIDs.append(i);
            }
        }
    }
    else
    {
        patchIDs = identity(patches.size());
    }
    patchIDs.shrink();

    labelHashSet calculatedPatches(patchIDs);

    label nTimeSteps = noiseParam.numberOfTimeSteps();
    instantList readTimes(nTimeSteps);


    scalar startTime = noiseParam.startTime();
    label startIndex = -1;
    forAll(allTimes, index)
    {
        if (allTimes[index].value() == startTime)
        {
            startIndex = index;
        }
    }

    if (startIndex == -1)
    {
        FatalErrorIn
        (
            "volumetricFFT::main(int argc, char *argv[])"
        )
            <<"Cannot find starting time: "<<startTime
            << exit(FatalError);
    }

    if
    (
        startIndex
        +(readTimes.size()-1)*noiseParam.timeRatio()
        >allTimes.size()
    )
    {

        FatalErrorIn
        (
            "volumetricFFT::main(int argc, char *argv[])"
        )
            << "first time-step, time ratio and number of samples "
            << "does no staisfy number of time-steps available: "<< allTimes.size()
            << exit(FatalError);
    }

    forAll(readTimes, timeI)
    {
        label index = startIndex + timeI*noiseParam.timeRatio();
        readTimes[timeI] = allTimes[index];
    }
    scalar endTime = readTimes[readTimes.size()-1].value();
    scalar deltaT = readTimes[1].value()-readTimes[0].value();

    scalar deltaF = 1.0/(endTime-startTime);

    scalar maxFreq = 1.0/(2.0*deltaT);
    scalar lowFreq = 20.0/(endTime-startTime);

    label ncells = mesh.nCells();

    if (returnReduce(ncells, sumOp<label>()))
    {
        scalar rho = 1.0;

        if (!noiseParam.compressible())
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties" ,
                    runTime.constant( ) ,
                    mesh ,
                    IOobject::MUST_READ ,
                    IOobject::NO_WRITE
                )
            ) ;
            rho = dimensionedScalar(transportProperties.lookup("rho")).value();
        }

        Info<<"Low frequency resolution (20 cycles) "<<lowFreq<<" Hz"<<endl;

        if (noiseParam.blockSize() > ncells)
        {
            noiseParam.blockSize() = ncells;
        }

        scalarField allTimeInternalData(0);
        if (noiseParam.internal())
        {
            allTimeInternalData.setSize(2*noiseParam.blockSize()*nTimeSteps, 0.0);
        }

        label nBoundaryFaces = 0 ;

        forAll(mesh.boundary( ), patchI)
        {
            if (!calculatedPatches.found(patchI))
            {
                continue;
            }

            nBoundaryFaces += mesh.boundary()[patchI].size();
        }

        scalarField allTimePatchData(2*nTimeSteps*nBoundaryFaces, 0.0);

        fftw_complex * fxSurf = reinterpret_cast <fftw_complex*>
            (allTimePatchData.begin());
        fftw_complex * fx = reinterpret_cast< fftw_complex*>
            (allTimeInternalData.begin());

        Info<< "Start time, endTime, number of time steps :"
             << startTime <<" "<<endTime<<" "<<nTimeSteps<< endl;

        label nBlocks = 1;
        if (noiseParam.blockSize() > 0)
        {
            nBlocks = max
            (
                1,
                ncells/noiseParam.blockSize()+((ncells%noiseParam.blockSize()) ? 1 : 0)
            );
        }

        DynamicList<scalar> octfm;
        DynamicList<DynamicList<scalar>> freqPerOneThirdOctave;
        if (noiseParam.calaculateOneThirdOctaveFFT())
        {
            calculateOctaveFrequencies
            (
                noiseParam,
                readTimes,
                deltaF,
                octfm,
                freqPerOneThirdOctave
            );
        }

        switchParallel sp(false);

        for
        (
            int blockI = 0;
            blockI < nBlocks;
            blockI++
        )
        {
            Pout<<"Calculating Block: " << blockI <<" of: "<<nBlocks<<endl;

            readPressure
            (
                mesh,
                noiseParam,
                readTimes,
                blockI,
                rho,
                calculatedPatches,
                allTimeInternalData,
                allTimePatchData,
                runTime
            );

            calculateDFT
            (
                noiseParam,
                blockI,
                nBoundaryFaces,
                fx,
                fxSurf
            );

            if (noiseParam.calaculateOneThirdOctaveFFT())
            {
                if (noiseParam.fUpper() > maxFreq)
                {
                    FatalErrorIn
                    (
                        "volumetricFFT::main(int argc, char *argv[])"
                    )
                        << "fUpper is higher than the sampled frequency: "
                        << maxFreq << " Hz"
                        << exit(FatalError);
                }

                if (noiseParam.fLower() == 0)
                {
                    FatalErrorIn
                    (
                        "volumetricFFT::main(int argc, char *argv[])"
                    )
                        << "fLower cannot be set to 0"
                        << exit(FatalError);
                }

                calculateOneThirds
                (
                   mesh,
                   noiseParam,
                   runTime,
                   fx,
                   fxSurf,
                   blockI,
                   ncells,
                   calculatedPatches,
                   readTimes,
                   deltaF,
                   octfm,
                   freqPerOneThirdOctave
                );
            }

            if (noiseParam.calculateEveryFrequency())
            {
                calculateAllFrequencies
                (
                    mesh,
                    noiseParam,
                    readTimes,
                    fx,
                    fxSurf,
                    blockI,
                    ncells,
                    patchIDs,
                    runTime
                );
            }

            if (noiseParam.calculateReverseFFT())
            {
                calculateReverseFFT
                (
                    mesh,
                    noiseParam,
                    readTimes,
                    blockI,
                    rho,
                    nBoundaryFaces,
                    ncells,
                    patchIDs,
                    deltaF,
                    fx,
                    fxSurf,
                    runTime
                );
            }
        }
    }
    Pout<<"Completed fft"<<endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
