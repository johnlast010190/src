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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/makeFvOption.H"
#include "FixedValueConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvOption(FixedValueConstraint, scalar);
makeFvOption(FixedValueConstraint, vector);
makeFvOption(FixedValueConstraint, sphericalTensor);
makeFvOption(FixedValueConstraint, symmTensor);
makeFvOption(FixedValueConstraint, tensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
void Foam::fv::FixedValueConstraint<Foam::vector>::constrain
(
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    vectorField values(cells_.size(), fieldValues_[fieldi]);
    for (direction cmptI = 0; cmptI < pTraits<vector>::nComponents; ++cmptI)
    {
        eqn.setValuesFromList(cmptI, cells_, values.component(cmptI));
    }
}


// ************************************************************************* //
