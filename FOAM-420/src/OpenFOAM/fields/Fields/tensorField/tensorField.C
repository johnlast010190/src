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
    (c) 2013 Esi Ltd.
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/Fields/tensorField/tensorField.H"
#include "fields/Fields/transformField/transformField.H"

#define TEMPLATE
#include "fields/Fields/Field/FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)

void inv(Field<tensor>& tf, const UList<tensor>& tf1)
{
    if (tf.empty())
    {
        return;
    }

    scalar scale = magSqr(tf1[0]);
    scalar smaller = SMALL*ROOTSMALL;

    Vector<bool> removeCmpts
    (
        magSqr(tf1[0].xx())/scale < smaller,
        magSqr(tf1[0].yy())/scale < smaller,
        magSqr(tf1[0].zz())/scale < smaller
    );

    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {

        if
        (
            (removeCmpts.x() && tf1[0].xx() != 0)
          ||(removeCmpts.y() && tf1[0].yy() != 0)
          ||(removeCmpts.z() && tf1[0].zz() != 0)
        )
        {
            WarningInFunction
                << "Employing stabilisation measures for tensor inversion." << nl
                << tab << "Accuracy of routines using this function are likely to "
                << "be affected." << nl << endl;
        }

        tensorField tf1Plus(tf1);

        if (removeCmpts.x())
        {
            tf1Plus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf1Plus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf1Plus += tensor(0,0,0,0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, hinv, tensor, tf1Plus)

        if (removeCmpts.x())
        {
            tf -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf -= tensor(0,0,0,0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, hinv, tensor, tf1)
    }
}

tmp<tensorField> inv(const UList<tensor>& tf)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result.ref(), tf);
    return result;
}

tmp<tensorField> inv(const tmp<tensorField>& tf)
{
    tmp<tensorField> tRes = New(tf);
    inv(tRes.ref(), tf());
    tf.clear();
    return tRes;
}

UNARY_FUNCTION(vector, tensor, eigenValues)
UNARY_FUNCTION(tensor, tensor, eigenVectors)

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)


template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const symmTensorField& stf
)
{
    tmp<tensorField> tRes(new tensorField(stf.size()));
    tensorField& res = tRes.ref();
    TFOR_ALL_F_OP_F(tensor, res, =, symmTensor, stf)
    return tRes;
}

template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<tensor>> ret = transformFieldMask<tensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/Fields/Field/undefFieldFunctionsM.H"

// ************************************************************************* //
