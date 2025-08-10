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
    (c) 2019 Esi Ltd.

Description


\*---------------------------------------------------------------------------*/

#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "faMatrices/faMatrix/faMatrix.H"
#include "finiteArea/d2dt2Schemes/faD2dt2Scheme/faD2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
d2dt2
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faD2dt2Scheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().d2dt2Scheme("d2dt2(" + vf.name() + ')')
    ).ref().famD2dt2(vf);
}


template<class Type>
tmp<faMatrix<Type>>
d2dt2
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faD2dt2Scheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().d2dt2Scheme
        (
            "d2dt2(" + rho.name() + ',' + vf.name() + ')'
        )
    ).ref().famD2dt2(rho, vf);
}


template<class Type>
tmp<faMatrix<Type>>
d2dt2
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faD2dt2Scheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().d2dt2Scheme
        (
            "d2dt2(" + rho.name() + ',' + vf.name() + ')'
        )
    ).ref().famD2dt2(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
