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
    (c) 2016, Esi Ltd

\*---------------------------------------------------------------------------*/

#include "noiseParameters.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::noiseParameters::noiseParameters
(
    const dictionary& dict
)
:
    blockSize_(dict.lookupOrDefault("blockSize", 10000000)),
    fLower_(dict.lookupOrDefault("fLower", 10.)),
    fUpper_(dict.lookupOrDefault("fUpper", 5000.)),
    hanning_(dict.lookupOrDefault("hanning", true)),
    hanningCoeff_(dict.lookupOrDefault("hanningCoeff", 1.)),
    calaculateOneThirdOctaveFFT_
    (dict.lookupOrDefault("calculateOneThirdOctaveFFT", false)),
    calculateEveryFrequency_
    (dict.lookupOrDefault("calculateEveryFrequency", false)),
    calculateReverseFFT_
    (dict.lookupOrDefault("calculateReverseFFT", false)),
    filterBands_
    (dict.lookupOrDefault("filterBands", vector2DField())),
    reverseFFTRatio_(dict.lookupOrDefault("reverseFFTRatio", 1)),
    startTime_(dict.lookupOrDefault("startTime", 0.)),
    numberOfTimeSteps_(readLabel(dict.lookup("numberOfTimeSteps", 1))),
    timeRatio_(dict.lookupOrDefault("timeRatio", 1)),
    internal_(dict.lookupOrDefault("internal", true)),
    compressible_(dict.lookupOrDefault("compressible", false)),
    binary_(dict.lookupOrDefault("binary", true)),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    patches_(dict.lookupOrDefault("patches", List<wordRe>()))
{
    Info<< nl
        << "Noise parameters" << nl
        << "----------------" << nl
        << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //

