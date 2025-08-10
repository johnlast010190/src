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

#include "finiteArea/fac/facDiv.H"
#include "faMesh/faMesh.H"
#include "finiteArea/fac/facEdgeIntegrate.H"
#include "finiteArea/divSchemes/faDivScheme/faDivScheme.H"
#include "finiteArea/convectionSchemes/faConvectionScheme/faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const GeometricField<Type, faePatchField, faEdgeMesh>& ssf
)
{
    return fac::edgeIntegrate(ssf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<GeometricField<Type, faePatchField, faEdgeMesh>>& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Div(fac::div(tssf()));
    tssf.clear();
    return Div;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::divScheme<Type>::New
    (
        vf.mesh(), vf.mesh().schemes().divScheme(name)
    ).ref().facDiv(vf);
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvvf,
    const word& name
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh>> Div
    (
        fac::div(tvvf(), name)
    );
    tvvf.clear();
    return Div;
}

template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::div(vf, "div("+vf.name()+')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvvf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh>> Div
    (
        fac::div(tvvf())
    );

    tvvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().schemes().divScheme(name)
    ).ref().facDiv(flux, vf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Div
    (
        fac::div(tflux(), vf, name)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Div
    (
        fac::div(flux, tvf(), name)
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Div
    (
        fac::div(tflux(), tvf(), name)
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::div
    (
        flux, vf, "div("+flux.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Div
    (
        fac::div(tflux(), vf)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Div
    (
        fac::div(flux, tvf())
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Div
    (
        fac::div(tflux(), tvf())
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
