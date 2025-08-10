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
    (c) 2015-2020 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "monolithicSolve.H"


template<>
Foam::solverPerformance Foam::monolithicSolve<Foam::scalar>::solve
(
    UPtrList<fvMatrix<scalar>>& matrices,
    const dictionary& solverControls,
    const bool meshMoved,
    const bool meshTopoChanged
)
{
    const label comm = Pstream::worldComm;

    UPtrList<GeometricField<scalar, fvPatchField, volMesh>> psi(matrices.size());
    forAll(matrices, i)
    {
        psi.set
        (
            i,
            &const_cast<GeometricField<scalar, fvPatchField, volMesh>&>
            (
                matrices[i].psi()
            )
        );
    }

    stringList fieldNames(psi.size());
    forAll(psi, i)
    {
        fieldNames[i] = psi[i].name();
    }
    word superName = makeSuperName(fieldNames);

    PtrList<scalarField> saveDiag(matrices.size());
    PtrList<scalarField> source(matrices.size());
    forAll(matrices, i)
    {
        matrices[i].boundaryManipulate(psi[i].boundaryFieldRef());

        saveDiag.set(i, new scalarField(matrices[i].diag()));
        source.set(i, new scalarField(matrices[i].source()));

        matrices[i].addBoundaryDiag(matrices[i].diag(), 0);
        matrices[i].addBoundarySource(source[i], matrices[i].psi(), false);
        // addBoundarySource will wrongly add the region coupled interfaces
        // because they do not report as coupled - this is compensated for
        // below
    }

    // Clear out relevant mesh objects associated to the super-mesh registry
    // if the mesh has changed since the last call
    Time& time = const_cast<Time&>(matrices[0].psi().mesh().time());
    if (meshMoved)
    {
        superMeshMoved(time);
    }
    if (meshTopoChanged)
    {
        superMeshTopoChanged(time);
    }

    // Concatenate all boundary fields as separate boundaries

    FieldField<Field, scalar> superBouCoeffs;
    FieldField<Field, scalar> superIntCoeffs;
    forAll(matrices, i)
    {
        const FieldField<Field, scalar>& bouCoeffs =
            matrices[i].boundaryCoeffs();
        const FieldField<Field, scalar>& intCoeffs =
            matrices[i].internalCoeffs();
        forAll(bouCoeffs, j)
        {
            const scalarField& pbc = bouCoeffs[j];
            const scalarField& pic = intCoeffs[j];
            superBouCoeffs.append(new Field<scalar>(pbc));
            superIntCoeffs.append(new Field<scalar>(pic));
            scalarField& psbc = superBouCoeffs.last();
            scalarField& psic = superIntCoeffs.last();

            fvPatchScalarField& pf = psi[i].boundaryFieldRef()[j];
            if (pf.regionCoupled())
            {
                regionCoupledFvPatchField<scalar>& rcpf =
                    refCast<regionCoupledFvPatchField<scalar>>
                    (
                        pf
                    );
                rcpf.regionCoupledBoundaryCoeffs
                (
                    matrices[i],
                    pbc,
                    pic,
                    psbc,
                    psic
                );
                scalarField fCorr( rcpf.faceCorr() );

                const labelUList& faceCells =
                    matrices[i].lduAddr().patchAddr(j);

                // Update the influence of the new int coeffs on the diagonal
                scalarField& diag( matrices[i].diag() );
                forAll(faceCells, l)
                {
                    diag[faceCells[l]] += (psic[l]-pic[l]);
                }

                // Update source for region-coupled correction
                forAll(faceCells, l)
                {
                    source[i][faceCells[l]] += pbc[l]*fCorr[l];
                }
            }
        }
    }

    // Make super matrix

    autoPtr<lduMesh> superMesh;
    autoPtr<lduMatrix> superMatrix;
    makeSuperMatrix(superMesh, superMatrix, matrices, time);

    // Concatenate all interface fields as separate boundaries

    PtrList<const lduInterfaceField> rawSuperInterfaceFields;
    UPtrList<const lduInterfaceField> superInterfaceFields;
    label k = 0;
    label offset = 0;
    forAll(matrices, i)
    {
        const lduInterfaceFieldPtrsList& interfaceFields =
            psi[i].boundaryField().scalarInterfaces();
        forAll(interfaceFields, j)
        {
            rawSuperInterfaceFields.resize(k+1);
            superInterfaceFields.resize(k+1);
            if (interfaceFields.set(j))
            {
                // See if the interface couples different regions
                label foreignOffset = offset; // offset of the other region mesh
                label foreignMatrix = i;
                if
                (
                    getCoupledMatrixAndOffset // Does nothing if not a coupled interface
                    (
                        interfaceFields[j].interface(),
                        matrices,
                        i,
                        j,
                        foreignMatrix,
                        foreignOffset
                    )
                )
                {
                    // Because region-coupled patches do not report as coupled,
                    // the boundary coefficient was incorrectly added to source
                    const labelUList& faceCells =
                        matrices[i].lduAddr().patchAddr(j);
                    const scalarField& pbc = matrices[i].boundaryCoeffs()[j];
                    const fvPatchScalarField& pf = psi[i].boundaryField()[j];
                    forAll(faceCells, l)
                    {
                        source[i][faceCells[l]] -= pbc[l]*pf[l];
                    }
                }
                rawSuperInterfaceFields.set
                (
                    k,
                    new offsetLduInterfaceField
                    (
                        interfaceFields[j],
                        interfaceFields[j].interface(),
                        superMesh->interfaces()[k],
                        offset,
                        matrices[i].lduAddr().size(),
                        foreignOffset,
                        (foreignMatrix == -1 ? 0 : matrices[foreignMatrix].lduAddr().size())
                    )
                );
                superInterfaceFields.set(k, &rawSuperInterfaceFields[k]);
            }
            k++;
        }
        offset += matrices[i].lduAddr().size(); // nCells
    }

    // Make super source and result fields
    scalarField superPsi(superMatrix->lduAddr().size());
    scalarField superSource(superMatrix->lduAddr().size());
    offset = 0;
    forAll(matrices, i)
    {
        forAll(psi[i], j)
        {
            superPsi[j+offset] = psi[i][j];
        }
        forAll(source[i], j)
        {
            superSource[j+offset] = source[i][j];
        }
        offset += psi[i].size();
    }

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        superName,
        superMatrix(),
        superBouCoeffs,
        superIntCoeffs,
        superInterfaceFields,
        solverControls
    )->solve(superPsi, superSource);

    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(comm));
    }

    // Unpack the result field
    offset = 0;
    forAll(psi, i)
    {
        psi[i].primitiveFieldRef() =
            scalarField::subField(superPsi, psi[i].size(), offset);
        offset += psi[i].size();
    }

    forAll(matrices, i)
    {
        matrices[i].diag() = saveDiag[i];
    }

    forAll(matrices, i)
    {
        psi[i].correctBoundaryConditions();
        psi[i].mesh().setSolverPerformance(psi[i].name(), solverPerf);
    }

    return solverPerf;
}
