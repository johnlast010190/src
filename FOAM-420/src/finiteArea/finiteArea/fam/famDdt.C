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


\*---------------------------------------------------------------------------*/

#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "faMatrices/faMatrix/faMatrix.H"
#include "finiteArea/ddtSchemes/faDdtScheme/faDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
ddt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme("ddt(" + vf.name() + ')')
    ).ref().famDdt(vf);
}


template<class Type>
tmp<faMatrix<Type>>
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme
        (
            "ddt(" + rho.name() + ',' + vf.name() + ')'
        )
    ).ref().famDdt(rho, vf);
}


template<class Type>
tmp<faMatrix<Type>>
ddt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddtScheme
        (
            "ddt(" + rho.name() + ',' + vf.name() + ')'
        )
    ).ref().famDdt(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
