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

\*---------------------------------------------------------------------------*/

#include "finiteArea/fac/facLnGrad.H"
#include "faMesh/faMesh.H"
#include "finiteArea/lnGradSchemes/lnGradScheme/lnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::lnGradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().lnGradScheme(name)
    )().lnGrad(vf);
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> LnGrad
    (
        fac::lnGrad(tvf(), name)
    );
    tvf.clear();
    return LnGrad;
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::lnGrad(vf, "lnGrad(" + vf.name() + ')');
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> LnGrad
    (
        fac::lnGrad(tvf())
    );
    tvf.clear();
    return LnGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
