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
    Abstract base class for edge interpolation schemes.

\*---------------------------------------------------------------------------*/

#include "interpolation/edgeInterpolation/edgeInterpolationScheme/edgeInterpolationScheme.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "fields/faPatchFields/faPatchField/faPatchFields.H"
#include "fields/faPatchFields/basic/coupled/coupledFaPatchField.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type>> edgeInterpolationScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::New(const faMesh&, Istream&)"
               " : constructing edgeInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTable_().find(schemeName);

    if (constructorIter == MeshConstructorTable_().end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << exit(FatalIOError);
    }

    return constructorIter->second(mesh, schemeData);
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type>> edgeInterpolationScheme<Type>::New
(
    const faMesh& mesh,
    const edgeScalarField& faceFlux,
    Istream& schemeData
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::New"
               "(const faMesh&, const edgeScalarField&, Istream&) : "
               "constructing edgeInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename MeshFluxConstructorTable::iterator constructorIter =
        MeshFluxConstructorTable_().find(schemeName);

    if (constructorIter == MeshFluxConstructorTable_().end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << exit(FatalIOError);
    }

    return constructorIter->second(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
edgeInterpolationScheme<Type>::~edgeInterpolationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the face-interpolate of the given cell field
//  with the given owner and neighbour weighting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas,
    const tmp<edgeScalarField>& tys
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::uncorrectedInterpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&, "
               "const tmp<edgeScalarField>&, "
               "const tmp<edgeScalarField>&) : "
               "interpolating areaTypeField from cells to faces "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();
    const edgeScalarField& ys = tys();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();
    const scalarField& y = ys.internalField();

    const faMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions(),
            faePatchField<Type>::calculatedType(),
            vf.boundaryField().patchTypes()
        )
    );
    GeometricField<Type, faePatchField, faEdgeMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    for (label fi=0; fi<P.size(); fi++)
    {
        // ZT, 22/Apr/2003
        const tensorField& curT = mesh.edgeTransformTensors()[fi];

        const tensor& Te = curT[0];
        const tensor& TP = curT[1];
        const tensor& TN = curT[2];

        sfi[fi] =
            transform
            (
                Te.T(),
                lambda[fi]*transform(TP, vfi[P[fi]])
              + y[fi]*transform(TN, vfi[N[fi]])
            );
    }


    // Interpolate across coupled patches using given lambdas and ys

    forAll(lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const faePatchScalarField& pY = ys.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            label size = vf.boundaryField()[pi].patch().size();
            label start = vf.boundaryField()[pi].patch().start();

            Field<Type> pOwnVf = vf.boundaryField()[pi].patchInternalField();
            Field<Type> pNgbVf = vf.boundaryField()[pi].patchNeighbourField();

            Field<Type>& pSf = sf.boundaryFieldRef()[pi];

            for (label i=0; i<size; i++)
            {
                const tensorField& curT =
                    mesh.edgeTransformTensors()[start + i];

                const tensor& Te = curT[0];
                const tensor& TP = curT[1];
                const tensor& TN = curT[2];

                pSf[i] =
                    transform
                    (
                        Te.T(),
                        pLambda[i]*transform(TP, pOwnVf[i])
                      + pY[i]*transform(TN, pNgbVf[i])
                    );
            }

//             sf.boundaryFieldRef()[pi] =
//                 pLambda*vf.boundaryField()[pi].patchInternalField()
//               + pY*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();
    tys.clear();

    return tsf;
}


//- Return the face-interpolate of the given cell field
//  with the given weigting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&, "
               "const tmp<edgeScalarField>&) : "
               "interpolating areaTypeField from cells to faces "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();

    const faMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions(),
            faePatchField<Type>::calculatedType(),
            vf.boundaryField().patchTypes()
        )
    );
    GeometricField<Type, faePatchField, faEdgeMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    for (label eI = 0; eI < P.size(); eI++)
    {
        // ZT, 22/Apr/2003
        const tensorField& curT = mesh.edgeTransformTensors()[eI];

        const tensor& Te = curT[0];
        const tensor& TP = curT[1];
        const tensor& TN = curT[2];

        sfi[eI] =
            transform
            (
                Te.T(),
                lambda[eI]*transform(TP, vfi[P[eI]])
              + (1 - lambda[eI])*transform(TN, vfi[N[eI]])
            );
    }


    // Interpolate across coupled patches using given lambdas

    forAll(lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            label size = vf.boundaryField()[pi].patch().size();
            label start = vf.boundaryField()[pi].patch().start();

            Field<Type> pOwnVf ( vf.boundaryField()[pi].patchInternalField() );
            Field<Type> pNgbVf ( vf.boundaryField()[pi].patchNeighbourField() );

            Field<Type>& pSf = sf.boundaryFieldRef()[pi];

            for (label i=0; i<size; i++)
            {
                const tensorField& curT =
                    mesh.edgeTransformTensors()[start + i];

                const tensor& Te = curT[0];
                const tensor& TP = curT[1];
                const tensor& TN = curT[2];

                pSf[i] =
                    transform
                    (
                        Te.T(),
                        pLambda[i]*transform(TP, pOwnVf[i])
                      + (1 - pLambda[i])*transform(TN, pNgbVf[i])
                    );
            }

//             tsf().boundaryFieldRef()[pi] =
//                 pLambda*vf.boundaryField()[pi].patchInternalField()
//              + (1 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


//- Return the euclidian edge-interpolate of the given area field
//  with the given weigting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
edgeInterpolationScheme<Type>::euclidianInterpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas
)
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::euclidianInterpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&, "
               "const tmp<edgeScalarField>&) : "
               "interpolating areaTypeField from cells to faces "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();

    const faMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions(),
            faePatchField<Type>::calculatedType(),
            vf.boundaryField().patchTypes()
        )
    );
    GeometricField<Type, faePatchField, faEdgeMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    for (label eI = 0; eI < P.size(); eI++)
    {
        sfi[eI] = lambda[eI]*vfi[P[eI]] + (1 - lambda[eI])*vfi[N[eI]];
    }


    // Interpolate across coupled patches using given lambdas

    forAll(lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            tsf.ref().boundaryFieldRef()[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
             + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


//- Return the face-interpolate of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&) : "
            << "interpolating areaTypeField from cells to faces"
            << endl;
    }

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf
        = interpolate(vf, weights(vf));

    if (corrected())
    {
        tsf.ref() += correction(vf);
    }

    return tsf;
}

//- Return the euclidian edge-interpolate of the given area field
//  without explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
edgeInterpolationScheme<Type>::euclidianInterpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    if (edgeInterpolation::debug)
    {
        Info<< "edgeInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, faPatchField, areaMesh>&) : "
            << "interpolating areaTypeField from cells to faces"
            << endl;
    }

    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf
        = euclidianInterpolate(vf, weights(vf));

    return tsf;
}

//- Return the face-interpolate of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
edgeInterpolationScheme<Type>::interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
) const
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tinterpVf
        = interpolate(tvf());
    tvf.clear();
    return tinterpVf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
