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
#include "finiteArea/laplacianSchemes/faLaplacianScheme/faLaplacianScheme.H"
#include "interpolation/edgeInterpolation/edgeInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type>
tmp<faMatrix<Type>>
laplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    edgeScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh().thisDb(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return fam::laplacian(Gamma, vf, name);
}


template<class Type>
tmp<faMatrix<Type>>
laplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    edgeScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh().thisDb(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return fam::laplacian
    (
        Gamma,
        vf,
        "laplacian(" + vf.name() + ')'
    );
}


template<class Type>
tmp<faMatrix<Type>>
laplacian
(
    const zero&,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return tmp<faMatrix<Type>>
    (
        new faMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<faMatrix<Type>>
laplacian
(
    const zero&,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<faMatrix<Type>>
    (
        new faMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<faMatrix<Type>>
laplacian
(
    const one&,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fam::laplacian(vf, name);
}


template<class Type>
tmp<faMatrix<Type>>
laplacian
(
    const one&,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::laplacian(vf);
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    const GeometricField<GType, faePatchField, faEdgeMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fam::laplacian(Gamma, vf, name);
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const GeometricField<GType, faePatchField, faEdgeMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fam::laplacian(Gamma, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<faMatrix<Type>>
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
    )().famLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const tmp<GeometricField<GType, faPatchField, areaMesh>>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type>> Laplacian(fam::laplacian(tgamma(), vf, name));
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const GeometricField<GType, faPatchField, areaMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const tmp<GeometricField<GType, faPatchField, areaMesh>>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> Laplacian(fam::laplacian(tgamma(), vf));
    tgamma.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<faMatrix<Type>>
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
    ).ref().famLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>& tgamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type>> tLaplacian = fam::laplacian(tgamma(), vf, name);
    tgamma.clear();
    return tLaplacian;
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const GeometricField<GType, faePatchField, faEdgeMesh>& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<faMatrix<Type>>
laplacian
(
    const tmp<GeometricField<GType, faePatchField, faEdgeMesh>>& tGamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam(fam::laplacian(tGamma(), vf));
    tGamma.clear();
    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
