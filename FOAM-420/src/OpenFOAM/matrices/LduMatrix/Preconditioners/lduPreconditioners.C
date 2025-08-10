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

#include "matrices/LduMatrix/Preconditioners/NoPreconditioner/NoPreconditioner.H"
#include "matrices/LduMatrix/Preconditioners/DiagonalPreconditioner/DiagonalPreconditioner.H"
#include "matrices/LduMatrix/Preconditioners/DILUPreconditioner/TDILUPreconditioner.H"
#include "fields/Fields/fieldTypes.H"

#define makeLduPreconditioners(Type, DType, LUType)                            \
                                                                               \
    makeLduPreconditioner(NoPreconditioner, Type, DType, LUType);              \
    makeLduSymPreconditioner(NoPreconditioner, Type, DType, LUType);           \
    makeLduAsymPreconditioner(NoPreconditioner, Type, DType, LUType);          \
                                                                               \
    makeLduPreconditioner(DiagonalPreconditioner, Type, DType, LUType);        \
    makeLduSymPreconditioner(DiagonalPreconditioner, Type, DType, LUType);     \
    makeLduAsymPreconditioner(DiagonalPreconditioner, Type, DType, LUType);    \
                                                                               \
    makeLduPreconditioner(TDILUPreconditioner, Type, DType, LUType);           \
    makeLduAsymPreconditioner(TDILUPreconditioner, Type, DType, LUType);

namespace Foam
{
    makeLduPreconditioners(scalar, scalar, scalar);
    makeLduPreconditioners(vector, scalar, scalar);
    makeLduPreconditioners(sphericalTensor, scalar, scalar);
    makeLduPreconditioners(symmTensor, scalar, scalar);
    makeLduPreconditioners(tensor, scalar, scalar);
};


// ************************************************************************* //
