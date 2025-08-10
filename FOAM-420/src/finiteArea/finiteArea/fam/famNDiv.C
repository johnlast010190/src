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
    (c) 2016-2017 Wikki Ltd

\*---------------------------------------------------------------------------*/

#include "finiteArea/fam/famNDiv.H"
#include "faMesh/faMesh.H"
#include "faMatrices/faMatrix/faMatrix.H"
#include "finiteArea/convectionSchemes/faConvectionScheme/faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const edgeScalarField& flux,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().schemes().divScheme(name)
    ).ref().famDiv(flux, vf);//TODO calculate normal
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type>> Div(fam::ndiv(tflux(), vf, name));
    tflux.clear();

    return Div;
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const edgeScalarField& flux,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::ndiv(flux, vf, "div("+flux.name()+','+vf.name()+')');
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> Div(fam::ndiv(tflux(), vf));
    tflux.clear();

    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
