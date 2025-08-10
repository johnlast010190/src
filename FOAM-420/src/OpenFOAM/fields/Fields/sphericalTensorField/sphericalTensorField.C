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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/Fields/sphericalTensorField/sphericalTensorField.H"
#include "fields/Fields/transformField/transformField.H"

#define TEMPLATE
#include "fields/Fields/Field/FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, inv)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)


template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const tensorField& tf
)
{
    return sph(tf);
}

template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<sphericalTensor>> ret =
        transformFieldMask<sphericalTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const symmTensorField& stf
)
{
    return sph(stf);
}

template<>
tmp<Field<sphericalTensor>> transformFieldMask<sphericalTensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<sphericalTensor>> ret =
        transformFieldMask<sphericalTensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/Fields/Field/undefFieldFunctionsM.H"

// ************************************************************************* //
