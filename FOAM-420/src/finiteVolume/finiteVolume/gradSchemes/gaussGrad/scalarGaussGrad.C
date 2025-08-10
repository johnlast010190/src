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

Description
    Specialisation of gaussGrad for scalars

\*---------------------------------------------------------------------------*/

#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, vector>> gaussGrad<scalar>::fvmGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();

    tmp<BlockLduSystem<vector, vector>> tbs
    (
        new BlockLduSystem<vector, vector>(mesh)
    );
    BlockLduSystem<vector, vector>& bs = tbs.ref();
    vectorField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    const surfaceVectorField mSf(fvc::applyFaceMask(mesh.Sf()));
    const vectorField& sf = mSf.internalField();

    tmp<surfaceScalarField> tweights = this->tinterpScheme_().weights(vf);
    const scalarField& w = tweights().internalField();

    l = -w*sf;
    u = l+sf;
    bs.negSumDiag();

    // Boundary contributions
    forAll(vf.boundaryField(), pI)
    {
        const fvPatchScalarField& pf = vf.boundaryField()[pI];
        const fvPatch& patch = pf.patch();
        const vectorField& pSf = mSf.boundaryField()[pI];
        const fvsPatchScalarField& pw = tweights().boundaryField()[pI];
        const labelList& fc = patch.faceCells();

        const scalarField iCoeffs(pf.valueInternalCoeffs(pw));
        const scalarField bCoeffs(pf.valueBoundaryCoeffs(pw));

        // Diag contribution
        forAll(pf, fI)
        {
            d[fc[fI]] += iCoeffs[fI]*pSf[fI];
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcUpper =
                bs.coupleUpper()[pI].asLinear();
            CoeffField<vector>::linearTypeField& pcLower =
                bs.coupleLower()[pI].asLinear();

            const vectorField pcl (iCoeffs*pSf);
            const vectorField pcu (bCoeffs*pSf);

            pcLower += pcl;
            pcUpper -= pcu;
        }
        else
        {
            forAll(pf, fI)
            {
                source[fc[fI]] -= bCoeffs[fI]*pSf[fI];
            }
        }
    }

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
