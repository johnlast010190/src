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

#include "finiteArea/ddtSchemes/steadyStateFaDdtScheme/steadyStateFaDdtScheme.H"
#include "finiteArea/fac/facDiv.H"
#include "faMatrices/faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensioned<Type> dt
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+dt.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensioned<Type> dt
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+dt.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt0("+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faePatchField, faEdgeMesh>& vf
)
{
    return tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "ddt0("+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt0("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "ddt0("+rho.name()+','+vf.name()+')',
                mesh()().time().timeName(),
                mesh()()
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                rho.dimensions()*vf.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );
}

template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );

    return tfam;
}

template<class Type>
tmp<edgeScalarField> steadyStateFaDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<edgeScalarField>
    (
        new edgeScalarField
        (
            IOobject
            (
                "meshPhi",
                mesh()().time().timeName(),
                mesh()(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("0", dimArea/dimTime, 0.0)
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
