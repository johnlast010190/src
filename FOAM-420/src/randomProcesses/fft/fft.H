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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::fft

Description
    Fast fourier transform using the fftw library.

    The complex transform field is returned in the field supplied.  The
    direction of transform is supplied as an argument (-1 = forward, 1 =
    reverse).  The dimensionality and organisation of the array of values
    in space is supplied in the nn indexing array.

Note
    The fftw library uses int only (no longs) for its dimensionality.

SourceFiles
    fft.C

\*---------------------------------------------------------------------------*/

#ifndef fft_H
#define fft_H

#include "fields/Fields/complexFields/complexFields.H"
#include "containers/Lists/UList/UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fft
{

public:

    enum transformDirection
    {
        FORWARD_TRANSFORM = -1,
        REVERSE_TRANSFORM = 1
    };


    static void transform
    (
        complexField& field,
        const UList<int>& nn,
        transformDirection fftDirection
    );


    static tmp<complexField> forwardTransform
    (
        const tmp<complexField>& field,
        const UList<int>& nn
    );


    static tmp<complexField> reverseTransform
    (
        const tmp<complexField>& field,
        const UList<int>& nn
    );


    static tmp<complexVectorField> forwardTransform
    (
        const tmp<complexVectorField>& field,
        const UList<int>& nn
    );


    static tmp<complexVectorField> reverseTransform
    (
        const tmp<complexVectorField>& field,
        const UList<int>& nn
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
