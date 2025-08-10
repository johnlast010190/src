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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
Su
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*su.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.source() -= mesh.S()*su.primitiveField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type>>
Su
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam = fam::Su(tsu(), vf);
    tsu.clear();
    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
Sp
(
    const areaScalarField& sp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.diag() += mesh.S()*sp.primitiveField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type>>
Sp
(
    const tmp<areaScalarField>& tsp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam = fam::Sp(tsp(), vf);
    tsp.clear();
    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
Sp
(
    const dimensionedScalar& sp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.diag() += mesh.S()*sp.value();

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
SuSp
(
    const areaScalarField& sp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.diag() += mesh.S()*max(sp.primitiveField(), scalar(0));

    fam.source() -= mesh.S()*min(sp.primitiveField(), scalar(0))
        *vf.primitiveField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type>>
SuSp
(
    const tmp<areaScalarField>& tsp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam = fam::SuSp(tsp(), vf);
    tsp.clear();
    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
