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
#include "lduRegisteredPrimitiveMesh/lduRegisteredPrimitiveMesh.H"
#include "mappedPatches/mappedLduInterface/mappedLduInterface.H"
#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"
#include "offsetLduInterface/offsetLduInterfaceField.H"
#include "meshes/MeshObject/MeshObject.H"


template<class Type>
bool Foam::monolithicSolve<Type>::getCoupledMatrixAndOffset
(
    const lduInterface& interface,
    const UPtrList<fvMatrix<Type>>& matrices,
    const label matrixI,
    const label patchI,
    label& foreignMatrix,
    label& foreignOffset
)
{
    // Don't couple through interfaces where the patch is coupled but not
    // the field - means this field is not coupled
    if
    (
        isA<mappedLduInterface>(interface)
     && matrices[matrixI].psi().boundaryField()[patchI].regionCoupled()
    )
    {
        const mappedLduInterface& rcInterface =
            refCast<const mappedLduInterface>
            (
                interface
            );
        foreignOffset = 0;
        bool found = false;
        forAll(matrices, j)
        {
            if (&matrices[j].psi().mesh() == &rcInterface.nbrMesh())
            {
                found = true;
                foreignMatrix = j;
                break;
            }
            foreignOffset += matrices[j].lduAddr().size();
        }
        if (!found)
        {
            WarningInFunction
                << "Coupled matrix not found for interface "
                << matrices[matrixI].psi().boundaryField()[patchI].patch().name()
                << endl;
            foreignOffset = 0;
            foreignMatrix = -1;
        }
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::monolithicSolve<Type>::makeSuperMatrix
(
    autoPtr<lduMesh>& superMesh,
    autoPtr<lduMatrix>& superMatrix,
    UPtrList<fvMatrix<Type>>& matrices,
    Time& time
)
{

    // Allocate lower, upper addressing

    label lSize = 0;
    label uSize = 0;
    label nCells = 0;
    forAll(matrices, i)
    {
        lSize += matrices[i].lduAddr().lowerAddr().size();
        uSize += matrices[i].lduAddr().upperAddr().size();
        nCells += matrices[i].lduAddr().size();
    }
    labelList sl(lSize);
    labelList su(uSize);

    // Concatenate lower, upper addressing

    label cellOffset = 0;
    label lowerOffset = 0;
    label upperOffset = 0;
    forAll(matrices, i)
    {
        // Append and adjust indices

        const labelList& l = matrices[i].lduAddr().lowerAddr();
        for (label j = 0; j < l.size(); j++)
        {
            sl[j+lowerOffset] = l[j]+cellOffset;
        }
        const labelList& u = matrices[i].lduAddr().upperAddr();
        for (label j = 0; j < u.size(); j++)
        {
            su[j+upperOffset] = u[j]+cellOffset;
        }

        // The losort addressing is created as needed, so we don't need to touch it.
        // Although we could actually create from the constituent matrices (if we
        // had access to the private members) since we aren't actually merging any
        // faces

        lowerOffset += l.size();
        upperOffset += u.size();
        cellOffset += matrices[i].lduAddr().size();
    }

    PtrList< const lduInterface > rawSuperInterfaces;
    lduInterfacePtrsList superInterfaces;
    lduSchedule superPatchSched;
    const label comm = Pstream::worldComm; // Make sure we select global comm, not mesh-specific

    // Concatenate all interfaces as separate boundaries
    label k = 0;
    cellOffset = 0;
    label lastPatchSchedNo = 0;
    forAll(matrices, i)
    {
        const lduInterfacePtrsList& interfaces =
            matrices[i].mesh().interfaces();
        label patchSchedOffset = lastPatchSchedNo;
        forAll(interfaces, j)
        {
            rawSuperInterfaces.resize(k+1);
            superInterfaces.resize(k+1);
            superPatchSched.resize(k+1);
            if (interfaces.set(j))
            {
                // See if the interface couples different regions
                label foreignOffset = cellOffset; // offset of the other region mesh
                label foreignMatrix = i;
                getCoupledMatrixAndOffset // Does nothing if not a coupled interface
                (
                    interfaces[j],
                    matrices,
                    i,
                    j,
                    foreignMatrix,
                    foreignOffset
                );

                // Adjust face-cell addressing
                rawSuperInterfaces.set
                (
                    k,
                    new offsetLduInterface
                    (
                        interfaces[j],
                        cellOffset,
                        matrices[i].lduAddr().size(),
                        foreignOffset,
                        (foreignMatrix == -1 ? 0 : matrices[foreignMatrix].lduAddr().size())
                    )
                );
                superInterfaces.set(k, &rawSuperInterfaces[k]);
            }
            superPatchSched[k] = matrices[i].lduAddr().patchSchedule()[j];
            superPatchSched[k].patch += patchSchedOffset;
            lastPatchSchedNo = superPatchSched[k].patch;
            k++;
        }
        cellOffset += matrices[i].lduAddr().size();
    }

    // Concatenate the region names for a unique name for this super-mesh
    // We do not want to keep the mesh in memory, but do keep a corresponding
    // objectRegistry for any mesh objects to be stored persistently
    string superName;
    for (auto mx : matrices)
    {
        superName += mx.psi().mesh().name();
    }
    if (!time.foundObject<objectRegistry>(superName))
    {
        time.store
        (
            new objectRegistry
            (
                IOobject
                (
                    superName,
                    time.timeName(),
                    time,
                    Foam::IOobject::NO_READ
                )
            )
        );
    }
    superMesh.set
    (
        new lduRegisteredPrimitiveMesh
        (
            nCells,
            sl,
            su,
            rawSuperInterfaces,
            superPatchSched,
            comm,
            time.lookupObject<objectRegistry>(superName)
        )
    );

    superMatrix.set
    (
        new lduMatrix(superMesh())
    );

    // Concatenate lower, upper and diag
    // Deal with any mix of diagonal, symmetric and asymmetric matrices

    bool hasLower = false;
    bool hasUpper = false;
    bool hasDiag = false;

    // Allocate super matrix to final size
    forAll(matrices, i)
    {
        hasLower = hasLower || matrices[i].hasLower();
        hasUpper = hasUpper || matrices[i].hasUpper();
        hasDiag = hasDiag || matrices[i].hasDiag();
    }
    if (hasLower)
    {
        superMatrix->lower();
    }
    if (hasUpper)
    {
        superMatrix->upper();
    }
    if (hasDiag)
    {
        superMatrix->diag();
    }

    // Populate

    label lOffset = 0;
    label uOffset = 0;
    label dOffset = 0;

    forAll(matrices, i)
    {

        label lSize = matrices[i].lduAddr().lowerAddr().size();
        label uSize = matrices[i].lduAddr().upperAddr().size();
        label dSize = matrices[i].lduAddr().size();
        if (hasLower)
        {
            scalarField& sl(superMatrix->lower());
            if (matrices[i].hasLower())
            {
                const scalarField& l(matrices[i].lower());
                forAll(l, j)
                {
                    sl[j+lOffset] = l[j];
                }
            }
            else if (matrices[i].hasUpper())
            {
                const scalarField& u(matrices[i].upper());
                forAll(u, j)
                {
                    sl[j+lOffset] = u[j];
                }
            }
            else
            {
                for (label j = 0; j < lSize; j++)
                {
                    sl[j+lOffset] = 0.0;
                }
            }
        }
        if (hasUpper)
        {
            scalarField& su(superMatrix->upper());
            if (matrices[i].hasUpper())
            {
                const scalarField& u(matrices[i].upper());
                forAll(u, j)
                {
                    su[j+uOffset] = u[j];
                }
            }
            else if (matrices[i].hasLower())
            {
                const scalarField& l(matrices[i].lower());
                forAll(l, j)
                {
                    su[j+uOffset] = l[j];
                }
            }
            else
            {
                for (label j = 0; j < uSize; j++)
                {
                    su[j+uOffset] = 0.0;
                }
            }
        }
        if (hasDiag)
        {
            scalarField& sd(superMatrix->diag());
            if (matrices[i].hasDiag())
            {
                const scalarField& d(matrices[i].diag());
                forAll(d, j)
                {
                    sd[j+dOffset] = d[j];
                }
            }
            else
            {
                for (label j = 0; j < dSize; j++)
                {
                    sd[j+dOffset] = 0.0;
                }
            }
        }

        lOffset += lSize;
        uOffset += uSize;
        dOffset += dSize;

    }

}


template <class Type>
void Foam::monolithicSolve<Type>::superMeshMoved(objectRegistry& obr)
{
    meshObject::clearUpto
    <
        lduMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(obr);
    meshObject::movePoints<lduMesh>(obr);
}


template <class Type>
void Foam::monolithicSolve<Type>::superMeshTopoChanged(objectRegistry& obr)
{
    meshObject::clearUpto
    <
        lduMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(obr);
    // Keep meshObjects that have an updateMesh
    // callback
    meshObject::clearUpto
    <
        lduMesh,
        TopologicalMeshObject,
        UpdateableMeshObject
    >
    (
        obr
    );
    //meshObject::updateMesh<lduMesh>(*this, mpm);
    // Too complicated to maintain and translate the map of changed cells to
    // the super mesh. Just remove and recreate the UpdateableMeshObjects
    meshObject::clear<lduMesh, UpdateableMeshObject>(obr);
}


//NOTE: There is a scalar specialisation of this function
template<class Type>
Foam::solverPerformance Foam::monolithicSolve<Type>::solve
(
    UPtrList<fvMatrix<Type>>& matrices,
    const dictionary& solverControls,
    const bool meshMoved,
    const bool meshTopoChanged
)
{
    const label comm = Pstream::worldComm;

    UPtrList<GeometricField<Type, fvPatchField, volMesh>> psi(matrices.size());
    forAll(matrices, i)
    {
        psi.set
        (
            i,
            &const_cast<GeometricField<Type, fvPatchField, volMesh>&>
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

    solverPerformance solverPerfVec
    (
        "solveMonolithic",
        superName
    );

    PtrList<scalarField> saveDiag(matrices.size());
    PtrList<Field<Type>> source(matrices.size());
    forAll(saveDiag, i)
    {
        matrices[i].boundaryManipulate(psi[i].boundaryFieldRef());

        saveDiag.set(i, new scalarField(matrices[i].diag()));
        source.set(i, new Field<Type>(matrices[i].source()));

        matrices[i].addBoundarySource(source[i], matrices[i].psi());
        // addBoundarySource will wrongly add the region coupled interfaces
        // because they do not report as coupled - this is compensated for
        // below
    }

    typename Type::labelType validComponents
    (
        pow
        (
            psi[0].mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );
    forAll(matrices, i)
    {
        if (psi[i].mesh().solutionD() != psi[0].mesh().solutionD())
        {
            FatalErrorInFunction
                << "Matrices have inconsistent solution directions."
                << endl
                << exit(FatalError);
        }
    }

    // Clear out relevant mesh objects associated to the super-mesh registry
    // if the mesh has changed since the last call
    Time& time =
        const_cast<Time&>(matrices[0].psi().mesh().time());
    if (meshMoved)
    {
        superMeshMoved(time);
    }
    if (meshTopoChanged)
    {
        superMeshTopoChanged(time);
    }

    // Convert region-coupled coefficients and get corrections

    FieldField<Field, Type> superBouCoeffs;
    FieldField<Field, Type> superIntCoeffs;
    PtrList<Field<Type>> superfCorr;
    forAll(matrices, i)
    {
        const FieldField<Field, Type>& bouCoeffs =
            matrices[i].boundaryCoeffs();
        const FieldField<Field, Type>& intCoeffs =
            matrices[i].internalCoeffs();
        forAll(bouCoeffs, j)
        {
            const Field<Type>& pbc = bouCoeffs[j];
            const Field<Type>& pic = intCoeffs[j];
            superBouCoeffs.append(new Field<Type>(pbc));
            superIntCoeffs.append(new Field<Type>(pic));
            Field<Type>& psbc = superBouCoeffs.last();
            Field<Type>& psic = superIntCoeffs.last();

            fvPatchField<Type>& pf = psi[i].boundaryFieldRef()[j];
            if (pf.regionCoupled())
            {
                regionCoupledFvPatchField<Type>& rcpf =
                    refCast<regionCoupledFvPatchField<Type>>
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
                superfCorr.set(j, rcpf.faceCorr().ptr());
            }
        }
    }


    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        PtrList<scalarField> psiCmpt(matrices.size());
        PtrList<scalarField> sourceCmpt(matrices.size());

        forAll(psiCmpt,i)
        {
            // copy field and source

            psiCmpt.set
            (
                i,
                new scalarField(psi[i].primitiveField().component(cmpt))
            );

            matrices[i].addBoundaryDiag(matrices[i].diag(), cmpt);

            sourceCmpt.set(i, new scalarField(source[i].component(cmpt)));
        }

        // Concatenate all boundary fields as separate boundaries

        FieldField<Field, scalar> superBouCoeffsCmpt;
        FieldField<Field, scalar> superIntCoeffsCmpt;
        forAll(matrices, i)
        {
            const FieldField<Field, Type>& bouCoeffs =
                matrices[i].boundaryCoeffs();
            const FieldField<Field, Type>& intCoeffs =
                matrices[i].internalCoeffs();
            forAll(bouCoeffs, j)
            {
                scalarField pbc(bouCoeffs[j].component(cmpt));
                scalarField pic(intCoeffs[j].component(cmpt));
                superBouCoeffsCmpt.append
                (
                    superBouCoeffs[j].component(cmpt).ptr()
                );
                superIntCoeffsCmpt.append
                (
                    superIntCoeffs[j].component(cmpt).ptr()
                );
                scalarField& psic = superIntCoeffsCmpt.last();

                fvPatchField<Type>& pf = psi[i].boundaryFieldRef()[j];
                if (pf.regionCoupled())
                {
                    scalarField fCorrCmpt(superfCorr[j].component(cmpt));

                    const labelUList& faceCells =
                        matrices[i].lduAddr().patchAddr(j);

                    // Update the influence of the new int coeffs on the diagonal
                    scalarField diagCmpt( matrices[i].diag().component(cmpt) );
                    forAll(faceCells, l)
                    {
                        diagCmpt[faceCells[l]] += (psic[l]-pic[l]);
                    }
                    matrices[i].diag().replace(cmpt, diagCmpt);

                    // Update source for region-coupled correction
                    forAll(faceCells, l)
                    {
                        sourceCmpt[i][faceCells[l]] += pbc[l]*fCorrCmpt[l];
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
                        const scalarField pbc
                        (
                            matrices[i].boundaryCoeffs()[j].component(cmpt)
                        );
                        const scalarField pf
                        (
                            psi[i].boundaryField()[j].component(cmpt)
                        );
                        forAll(faceCells, l)
                        {
                            sourceCmpt[i][faceCells[l]] -= pbc[l]*pf[l];
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
                            (
                                foreignMatrix == -1 ?
                                0 :
                                matrices[foreignMatrix].lduAddr().size()
                            )
                        )
                    );
                    superInterfaceFields.set(k, &rawSuperInterfaceFields.last());
                }
                k++;
            }
            offset += matrices[i].lduAddr().size(); // nCells
        }

        // Make super source and result fields

        scalarField superPsiCmpt(superMatrix->lduAddr().size());
        scalarField superSourceCmpt(superMatrix->lduAddr().size());
        offset = 0;
        forAll(matrices, i)
        {
            forAll(psiCmpt[i], j)
            {
                superPsiCmpt[j+offset] = psiCmpt[i][j];
            }
            forAll(sourceCmpt[i], j)
            {
                superSourceCmpt[j+offset] = sourceCmpt[i][j];
            }
            offset += psi[i].size();
        }

        // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
        // bouCoeffsCmpt for the explicit part of the coupled boundary
        // conditions
        superMatrix->initMatrixInterfaces
        (
            true,
            superBouCoeffsCmpt,
            superInterfaceFields,
            superPsiCmpt,
            superSourceCmpt,
            cmpt
        );

        superMatrix->updateMatrixInterfaces
        (
            true,
            superBouCoeffsCmpt,
            superInterfaceFields,
            superPsiCmpt,
            superSourceCmpt,
            cmpt
        );


        solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            superName + pTraits<Type>::componentNames[cmpt],
            superMatrix(),
            superBouCoeffsCmpt,
            superIntCoeffsCmpt,
            superInterfaceFields,
            solverControls
        )->solve(superPsiCmpt, superSourceCmpt, cmpt);

        if (solverPerformance::debug)
        {
            solverPerf.print(Info.masterStream(comm));
        }

        solverPerfVec = max(solverPerfVec, solverPerf);
        solverPerfVec.solverName() = solverPerf.solverName();

        // Unpack the result field
        offset = 0;
        forAll(psiCmpt, i)
        {
            psi[i].primitiveFieldRef().replace
            (
                cmpt,
                scalarField::subList(superPsiCmpt, psiCmpt[i].size(), offset)
            );
            offset += psiCmpt[i].size();
        }

        forAll(saveDiag, i)
        {
            matrices[i].diag() = saveDiag[i];
        }
    }

    forAll(psi, i)
    {
        psi[i].correctBoundaryConditions();
        psi[i].mesh().setSolverPerformance(psi[i].name(), solverPerfVec);
    }

    return solverPerfVec;
}


template<class Type>
Foam::word Foam::monolithicSolve<Type>::makeSuperName
(
    const stringList& names
)
{
    // Concatenate the unique names of constituent fields
    stringList sortedNames(names);

    // Sort alphabetically to get a predictable ordering
    sort(sortedNames);

    DynamicList<string> uniqueNames;
    uniqueNames.append(sortedNames[0]);
    for (label i = 1; i < sortedNames.size(); i++)
    {
        if (sortedNames[i] != sortedNames[i-1])
        {
            uniqueNames.append(sortedNames[i]);
        }
    }

    string superName;
    forAll(uniqueNames, i)
    {
        superName += uniqueNames[i];
        if (i < uniqueNames.size()-1)
        {
            superName += "_";
        }
    }
    return superName;
}


template<class Type>
Foam::solverPerformance Foam::monolithicSolve<Type>::solve
(
    UPtrList<fvMatrix<Type>>& matrices,
    const solution& solutionControls,
    const bool finalIter,
    const bool meshMoved,
    const bool meshTopoChanged
)
{
    stringList names(matrices.size());
    forAll(matrices, i)
    {
        names[i] = matrices[i].psi().name();
    }
    word superName = makeSuperName(names);
    const dictionary& solverControls =
        solutionControls.solverDict(superName + (finalIter ? "Final" : ""));
    return
        solve(matrices, superName, solverControls, meshMoved, meshTopoChanged);
}
