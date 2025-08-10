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
    Simple central-difference lnGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "finiteArea/lnGradSchemes/correctedLnGrad/correctedLnGrad.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "interpolation/edgeInterpolation/schemes/linear/linearEdgeInterpolation.H"
#include "finiteArea/fac/facGrad.H"
#include "finiteArea/gradSchemes/gaussFaGrad/gaussFaGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fa::correctedLnGrad<Type>::~correctedLnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::faEdgeMesh>>
Foam::fa::correctedLnGrad<Type>::fullGradCorrection
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    const faMesh& mesh = this->mesh();

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tssf =
        mesh.nonOrthCorrectionVectors()
      & linearEdgeInterpolation<typename outerProduct<vector, Type>::type>(mesh).interpolate
        (
            gradScheme<Type>::New
            (
                mesh,
                mesh.schemes().gradScheme("grad(" + vf.name() + ')')
            )().grad(vf)
        );
    tssf.ref().rename("lnGradCorr(" + vf.name() + ')');

    return tssf;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::faEdgeMesh>>
Foam::fa::correctedLnGrad<Type>::correction
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    const faMesh& mesh = this->mesh();

    // construct GeometricField<Type, faePatchField, faEdgeMesh>
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tssf
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "lnGradCorr("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, faePatchField, faEdgeMesh>& ssf = tssf.ref();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        ssf.replace
        (
            cmpt,
            correctedLnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
