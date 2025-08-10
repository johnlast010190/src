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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2018-2023 ESI

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Constant/Constant.H"
#include "primitives/functions/Function1/Zero/ZeroConstant.H"
#include "primitives/functions/Function1/One/OneConstant.H"
#include "primitives/functions/Function1/PolynomialEntry/PolynomialEntry.H"
#include "primitives/functions/Function1/NSRDSfunc0/NSRDSfunc0.H"
#include "primitives/functions/Function1/NSRDSfunc1/NSRDSfunc1.H"
#include "primitives/functions/Function1/NSRDSfunc2/NSRDSfunc2.H"
#include "primitives/functions/Function1/NSRDSfunc3/NSRDSfunc3.H"
#include "primitives/functions/Function1/NSRDSfunc4/NSRDSfunc4.H"
#include "primitives/functions/Function1/NSRDSfunc5/NSRDSfunc5.H"
#include "primitives/functions/Function1/NSRDSfunc6/NSRDSfunc6.H"
#include "primitives/functions/Function1/NSRDSfunc7/NSRDSfunc7.H"
#include "primitives/functions/Function1/NSRDSfunc14/NSRDSfunc14.H"
#include "primitives/functions/Function1/Sine/Sine.H"
#include "primitives/functions/Function1/Square/Square.H"
#include "primitives/functions/Function1/CSV/CSV.H"
#include "primitives/functions/Function1/Table/Table.H"
#include "primitives/functions/Function1/TableFile/TableFile.H"
#include "primitives/functions/Function1/Scale/Scale.H"
#include "primitives/functions/Function1/SelectEntry/SelectEntry.H"
#include "primitives/functions/Function1/xyzPolynomial/xyzPolynomial.H"

#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction1s(Type)                                                   \
    makeFunction1(Type);                                                       \
    makeFunction1Type(Constant, Type);                                         \
    makeFunction1Type(ZeroConstant, Type);                                     \
    makeFunction1Type(OneConstant, Type);                                      \
    makeFunction1Type(Polynomial, Type);                                       \
    makeFunction1Type(xyzPolynomial, Type);                                    \
    makeFunction1Type(Sine, Type);                                             \
    makeFunction1Type(Square, Type);                                           \
    makeFunction1Type(CSV, Type);                                              \
    makeFunction1Type(Table, Type);                                            \
    makeFunction1Type(TableFile, Type);                                        \
    makeFunction1Type(Scale, Type);                                            \
    makeFunction1Type(Select, Type);                                           \
    makeFunction1Type(NSRDSfunc0, Type);                                       \
    makeFunction1Type(NSRDSfunc1, Type);                                       \
    makeFunction1Type(NSRDSfunc2, Type);                                       \
    makeFunction1Type(NSRDSfunc3, Type);                                       \
    makeFunction1Type(NSRDSfunc4, Type);                                       \
    makeFunction1Type(NSRDSfunc5, Type);                                       \
    makeFunction1Type(NSRDSfunc6, Type);                                       \
    makeFunction1Type(NSRDSfunc7, Type);                                       \
    makeFunction1Type(NSRDSfunc14, Type);



namespace Foam
{
    makeFunction1(label);
    makeFunction1Type(Constant, label);
    // Polynomial functions and interpolation do evaluate to label
    // Instead evaluate a scalar and convert to label as appropriate

    makeFunction1s(scalar);
    makeFunction1s(vector);
    makeFunction1s(sphericalTensor);
    makeFunction1s(symmTensor);
    makeFunction1s(tensor);

    makeFunction1(bool);
    makeFunction1Type(Table, bool);
    makeFunction1Type(Constant, bool);
}


// ************************************************************************* //
