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
    Abstract base class for lnGrad schemes.

\*---------------------------------------------------------------------------*/

#include "finiteArea/fa/fa.H"
#include "finiteArea/lnGradSchemes/lnGradScheme/lnGradScheme.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "containers/HashTables/HashTable/HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<lnGradScheme<Type>> lnGradScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (fa::debug)
    {
        Info<< "lnGradScheme<Type>::New(const faMesh&, Istream&)"
               " : constructing lnGradScheme<Type>"
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
lnGradScheme<Type>::~lnGradScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the face-lnGrad of the given cell field
//  with the given weigting factors
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGradScheme<Type>::lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tdeltaCoeffs,
    const word& lnGradName
)
{
    const faMesh& mesh = vf.mesh();

    // construct GeometricField<Type, faePatchField, faEdgeMesh>
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tssf
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                lnGradName + vf.name() + ')',
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*tdeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, faePatchField, faEdgeMesh>& ssf = tssf.ref();

    // set reference to difference factors array
    const scalarField& deltaCoeffs = tdeltaCoeffs().internalField();

    // owner/neighbour addressing
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        ssf[faceI] =
            deltaCoeffs[faceI]*(vf[neighbour[faceI]] - vf[owner[faceI]]);
    }

    forAll(vf.boundaryField(), patchI)
    {
        const faPatchField<Type>& pvf = vf.boundaryField()[patchI];

            ssf.boundaryFieldRef()[patchI] = pvf.snGrad();
    }

    return tssf;
}


//- Return the face-lnGrad of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGradScheme<Type>::lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf
        = lnGrad(vf, deltaCoeffs(vf));

    if (corrected())
    {
        tsf.ref() += correction(vf);
    }

    return tsf;
}


//- Return the face-lnGrad of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
lnGradScheme<Type>::lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
) const
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf
    (
        lnGrad(tvf())
    );

    tsf.clear();
    return tsf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
