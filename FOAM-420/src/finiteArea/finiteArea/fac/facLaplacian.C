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
    Namespace of functions to calculate explicit derivatives.
    Time derivatives are calculated using Euler-implicit, backward differencing
    or Crank-Nicholson. Spatial derivatives are calculated using Gauss' Theorem.

\*---------------------------------------------------------------------------*/

#include "finiteArea/fac/facLaplacian.H"
#include "faMesh/faMesh.H"
#include "finiteArea/laplacianSchemes/faLaplacianScheme/faLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type, scalar>::New
    (
        vf.mesh(),
        vf.mesh().schemes().laplacianScheme(name)
    ).ref().facLaplacian(vf);
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian(vf, "laplacian(" + vf.name() + ')');
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tvf())
    );
    tvf.clear();
    return Laplacian;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    GeometricField<GType, faePatchField, faEdgeMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fac::laplacian(Gamma, vf, name);
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<GType, faePatchField, faEdgeMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fac::laplacian(Gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faPatchField, areaMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().schemes().laplacianScheme(name)
    ).ref().facLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<GType, faPatchField, areaMesh>>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faPatchField, areaMesh>& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<GType, faPatchField, areaMesh>>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faPatchField, areaMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<GType, faPatchField, areaMesh>>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian
    (
        tgamma,
        vf,
        "laplacian(" + tgamma().name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faPatchField, areaMesh>& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    return fac::laplacian
    (
        gamma,
        tvf,
        "laplacian(" + gamma.name() + ',' + tvf().name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<GType, faPatchField, areaMesh>>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    return fac::laplacian
    (
        tgamma,
        tvf,
        "laplacian(" + tgamma().name() + ',' + tvf().name() + ')'
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faePatchField, faEdgeMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().schemes().laplacianScheme(name)
    ).ref().facLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faePatchField, faEdgeMesh>& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>> laplacian
(
    const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faePatchField, faEdgeMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tgamma(), vf)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>>
laplacian
(
    const GeometricField<GType, faePatchField, faEdgeMesh>& gamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, faPatchField, areaMesh>> laplacian
(
    const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>& tgamma,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> Laplacian
    (
        fac::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
