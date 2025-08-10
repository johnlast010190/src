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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Application
    noise

Group
    grpPostProcessingUtilities

Description
    Utility to perform noise analysis of pressure data.

    The utility provides a light wrapper around the run-time selectable
    noise model.  Current options include:
    - point, and
    - surface noise.

Usage
    \verbatim
    noiseModel      surfaceNoise; // pointNoise

    surfaceNoiseCoeffs
    {
        windowModel     Hanning;

        HanningCoeffs
        {
            // Window overlap percentage
            overlapPercent  50;
            symmetric       yes;
            extended        yes;

            // Optional number of windows, default = all available
            nWindow         5;
        }


        // Input file
        inputFile   "postProcessing/faceSource1/surface/patch/patch.case";

        // Surface reader
        reader      ensight;

        // Surface writer
        writer      ensight;

        // Collate times for ensight output
        // - ensures geometry is only written once
        writeOptions
        {
            ensight
            {
                collateTimes true;
            }
        }

        // Number of samples in sampling window
        // Must be a power of 2, default = 2^16 (=65536)
        N               4096;

        // Write interval for FFT data, default = 1
        fftWriteInterval 100;
    }
    \endverbatim


See also
    noiseFFT.H
    noiseModel.H
    windowModel.H

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "noise/noiseModels/noiseModel/noiseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    #include "include/addDictOption.H"
    #include "include/setRootCase.H"
    #include "include/createTime.H"

    const word dictName("noiseDict");
    #include "include/setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << "\n" << endl;

    IOdictionary dict(dictIO);

    autoPtr<noiseModel> model(noiseModel::New(dict));
    model->calculate();

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
