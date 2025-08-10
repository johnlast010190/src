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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Description
    Simple central-difference lnGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "finiteArea/lnGradSchemes/correctedLnGrad/correctedLnGrad.H"
#include "faMesh/faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeLnGradScheme(correctedLnGrad)

namespace Foam
{
namespace fa
{
    template<>
    tmp<Foam::edgeScalarField>
    correctedLnGrad<scalar>::correction
    (
        const areaScalarField& vsf
    ) const
    {
        return fullGradCorrection(vsf);
    }


    template<>
    tmp<edgeVectorField>
    correctedLnGrad<vector>::correction
    (
        const areaVectorField& vvf
    ) const
    {
        return fullGradCorrection(vvf);
    }
}
}

// ************************************************************************* //
