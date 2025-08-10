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

Class
    Foam::noiseParameters

Description
    Simple container for keeping all noise specific information

SourceFiles
    noiseParameters.C
\*---------------------------------------------------------------------------*/


#ifndef noiseParameters_H
#define noiseParameters_H

#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "fields/Fields/vector2DField/vector2DField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class noiseParameters Declaration
\*---------------------------------------------------------------------------*/

class noiseParameters
{

private:

    label blockSize_;

    scalar fLower_;

    scalar fUpper_;

    Switch hanning_;

    scalar hanningCoeff_;

    Switch calaculateOneThirdOctaveFFT_;

    Switch calculateEveryFrequency_;

    Switch calculateReverseFFT_;

    vector2DField filterBands_;

    label reverseFFTRatio_;

    scalar startTime_;

    label numberOfTimeSteps_;

    label timeRatio_;

    Switch internal_;

    Switch compressible_;

    Switch binary_;

    word pName_;

    List<wordRe> patches_;

public:

    // Constructors

    //- Construct from dictionary
    noiseParameters(const dictionary& dict);

    // Member Functions

    //Split internal cells in blocks of given size to reduce memory usage
    label blockSize() const
    {
        return blockSize_;
    }

    label& blockSize()
    {
        return blockSize_;
    }

    //Lower frequency for calculating octave bands
    scalar fLower() const
    {
        return fLower_;
    }

    //Upper frequency for calculating octave bands
    scalar fUpper() const
    {
        return fUpper_;
    }

    //Whether of not to apply a Hanning function
    Switch hanning() const
    {
        return hanning_;
    }

    //Set Hanning coefficient
    scalar hanningCoeff() const
    {
        return hanningCoeff_;
    }

    //Whether to output third-octaves
    Switch calaculateOneThirdOctaveFFT() const
    {
        return calaculateOneThirdOctaveFFT_;
    }

    //Whether to output for all frequencies
    Switch calculateEveryFrequency() const
    {
        return calculateEveryFrequency_;
    }

    //Whether to output reverse fourier transforms
    Switch calculateReverseFFT() const
    {
        return calculateReverseFFT_;
    }

    //Set filter bands for reverse transforms
    vector2DField filterBands() const
    {
        return filterBands_;
    }

    //Ratio of FFT data for inverse transforms
    label reverseFFTRatio() const
    {
        return reverseFFTRatio_;
    }


    //starting time
    scalar startTime() const
    {
        return startTime_;
    }

    //Number of time steps to parse
    label numberOfTimeSteps() const
    {
        return numberOfTimeSteps_;
    }

    //ratio used for time steps
    label timeRatio() const
    {
        return timeRatio_;
    }

    //Whether to calculate for interior cells
    Switch internal() const
    {
        return internal_;
    }

    //Whether case is compressible
    Switch compressible() const
    {
        return compressible_;
    }

    //Whether to write in binary
    Switch binary() const
    {
        return binary_;
    }

    //Name of pressure field
    word pName() const
    {
        return pName_;
    }

    //list of patches to calculate over
    List<wordRe> patches() const
    {
        return patches_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
