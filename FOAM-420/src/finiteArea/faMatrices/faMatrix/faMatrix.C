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
    Finite-Area matrix

\*---------------------------------------------------------------------------*/
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "fields/faPatchFields/basic/calculated/calculatedFaPatchFields.H"
#include "fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchFields.H"
#include "fields/faPatchFields/basic/coupled/coupledFaPatchFields.H"
#include "containers/Lists/UIndirectList/UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void faMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    FOAM_ASSERT(addr.size() == pf.size()) {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] += pf[faceI];
    }
}


template<class Type>
template<class Type2>
void faMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void faMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] -= pf[faceI];
    }
}


template<class Type>
template<class Type2>
void faMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
void faMatrix<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            internalCoeffs_[patchI].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void faMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            cmptAv(internalCoeffs_[patchI]),
            diag
        );
    }
}


template<class Type>
void faMatrix<Type>::addBoundarySource
(
    Field<Type>& source,
    const GeometricField<Type, faPatchField, areaMesh>& psi,
    const bool couples
) const
{
    forAll(psi.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi.boundaryField()[patchI];
        const Field<Type>& pbc = boundaryCoeffs_[patchI];

        if (!ptf.coupled())
        {
            addToInternalField(lduAddr().patchAddr(patchI), pbc, source);
        }
        else if (couples)
        {
            tmp<Field<Type>> tpnf = ptf.patchNeighbourField();
            const Field<Type>& pnf = tpnf();

            const labelUList& addr = lduAddr().patchAddr(patchI);

            forAll(addr, facei)
            {
                source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
            }
        }
    }
}

template<class Type>
template<template<class> class ListType>
void faMatrix<Type>::setValuesFromList
(
        const labelUList& faceLabels,
        const ListType<Type>& values
)
{
const faMesh& mesh = psi_.mesh();

const labelListList& edges = mesh.patch().faceEdges();
const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

scalarField& Diag = diag();
Field<Type>& psi =
    const_cast
    <
        GeometricField<Type, faPatchField, areaMesh>&
    >(psi_).primitiveFieldRef();

forAll(faceLabels, i)
{
    const label facei = faceLabels[i];
    const Type& value = values[i];

    psi[facei] = value;
    source_[facei] = value*Diag[facei];
    if (symmetric() || asymmetric())
    {
        const labelList& f= edges[facei];
        label nBound = 0;
        forAll(f, j)
        {
            const label edgei = f[j];

            if (mesh.isInternalEdge(edgei))
            {
                if (symmetric())
                {
                    if (facei == own[edgei])
                    {
                        source_[nei[edgei]] -= upper()[edgei]*value;
                    }
                    else
                    {
                        source_[own[edgei]] -= upper()[edgei]*value;
                    }

                    upper()[edgei] = 0.0;
                }
                else
                {
                    if (facei == own[edgei])
                    {
                        source_[nei[edgei]] -= lower()[edgei]*value;
                    }
                    else
                    {
                        source_[own[edgei]] -= upper()[edgei]*value;
                    }

                    upper()[edgei] = 0.0;
                    lower()[edgei] = 0.0;
                }
            }
            else
            {
                label patchi = mesh.boundary().whichPatch(edgei);
                if (internalCoeffs_[patchi].size())
                {
                    label patchFacei =
                        mesh.boundary()[patchi].whichEdge(edgei);

                    internalCoeffs_[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    boundaryCoeffs_[patchi][patchFacei] =
                        pTraits<Type>::zero;
                }
                nBound++;
            }
        }
        //reset source and diagonal if all cell faces boundary
        if (nBound == f.size())
        {
            Diag[facei] = 1.0;
            source_[facei] = value*Diag[facei];
        }
    }
}
}




// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
faMatrix<Type>::faMatrix
(
    const GeometricField<Type, faPatchField, areaMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), pTraits<Type>::zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(nullptr),
    solvingComponent(0)
{
    if (debug)
    {
        Info<< "faMatrix<Type>(GeometricField<Type, faPatchField, areaMesh>&,"
               " const dimensionSet&) : "
               "constructing faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

    // Update the boundary coefficients of psi without changing its event No.
    GeometricField<Type, faPatchField, areaMesh>& psiRef =
       const_cast<GeometricField<Type, faPatchField, areaMesh>&>(psi_);

    label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryFieldRef().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
faMatrix<Type>::faMatrix(const faMatrix<Type>& fam)
:
    refCount(),
    lduMatrix(fam),
    psi_(fam.psi_),
    dimensions_(fam.dimensions_),
    source_(fam.source_),
    internalCoeffs_(fam.internalCoeffs_),
    boundaryCoeffs_(fam.boundaryCoeffs_),
    faceFluxCorrectionPtr_(nullptr),
    solvingComponent(fam.solvingComponent)
{
    if (debug)
    {
        Info<< "faMatrix<Type>::faMatrix(const faMatrix<Type>&) : "
            << "copying faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (fam.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, faePatchField, faEdgeMesh>
        (
            *(fam.faceFluxCorrectionPtr_)
        );
    }
}


template<class Type>
faMatrix<Type>::faMatrix
(
    const GeometricField<Type, faPatchField, areaMesh>& psi,
    Istream& is
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(is),
    source_(is),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(nullptr),
    solvingComponent(0)
{
    if (debug)
    {
        Info<< "faMatrix<Type>(GeometricField<Type, faPatchField, areaMesh>&,"
               " Istream&) : "
               "constructing faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

}


template<class Type>
faMatrix<Type>::~faMatrix()
{
    if (debug)
    {
        Info<< "faMatrix<Type>::~faMatrix<Type>() : "
            << "destroying faMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (faceFluxCorrectionPtr_)
    {
        delete faceFluxCorrectionPtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// turn off solution of given cells, retaining current value
// and eliminating corresponding equation from the matrix
template<class Type>
void faMatrix<Type>::fixValues
(
    const labelList& faceLabels
)
{
    const faMesh& mesh = psi_.mesh();

    const labelListList& edges = mesh.patch().faceEdges();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    scalarField& Diag = diag();

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        source_[facei] = psi_[facei]*Diag[facei];

        if (symmetric() || asymmetric())
        {
            const labelList& f= edges[facei];

            forAll(f, j)
            {
                label edgei = f[j];

                if (mesh.isInternalEdge(edgei))
                {
                    if (symmetric())
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= upper()[edgei]*psi_[facei];
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*psi_[facei];
                        }

                        upper()[edgei] = 0.0;
                    }
                    else
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= lower()[edgei]*psi_[facei];
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*psi_[facei];
                        }

                        upper()[edgei] = 0.0;
                        lower()[edgei] = 0.0;
                    }
                }
                else
                {
                    label patchi = mesh.boundary().whichPatch(edgei);

                    if (internalCoeffs_[patchi].size())
                    {
                        label patchEdgei =
                            mesh.boundary()[patchi].whichEdge(edgei);

                        internalCoeffs_[patchi][patchEdgei] =
                            pTraits<Type>::zero;

                        boundaryCoeffs_[patchi][patchEdgei] =
                            pTraits<Type>::zero;
                    }
                }
            }
        }
    }
}

template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const UList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}

template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& cellLabels,
    const UIndirectList<Type>& values
)
{
    this->setValuesFromList(cellLabels, values);
}

// Set reference level for solution
template<class Type>
void faMatrix<Type>::setReference
(
    const label facei,
    const Type& value,
    const bool forceReference
)
{
    if ((forceReference || psi_.needReference()) && facei >= 0)
    {
        source()[facei] += diag()[facei]*value;
        diag()[facei] += diag()[facei];
    }

}

template<class Type>
void faMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    if (debug)
    {
        InfoIn("faMatrix<Type>::relax(const scalar alpha)")
            << "Relaxing " << psi_.name() << " by " << alpha
            << endl;
    }

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), 0.0);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                const Field<Type>& pCoeffs = boundaryCoeffs_[patchI];

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                forAll(pa, face)
                {
                    D[pa[face]] += component(iCoeffs[face], 0);
                    sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
                }
            }
            else
            {
                // For non-coupled boundaries add the maximum magnitude diagonal
                // contribution to ensure stability
                forAll(pa, face)
                {
                    D[pa[face]] += cmptMax(cmptMag(iCoeffs[face]));
                }
            }
        }
    }

    if (debug)
    {
        // Calculate amount of non-dominance.
        label nNon = 0;
        scalar maxNon = 0.0;
        scalar sumNon = 0.0;
        forAll(D, facei)
        {
            scalar d = (sumOff[facei] - D[facei])/mag(D[facei]);

            if (d > 0)
            {
                nNon++;
                maxNon = max(maxNon, d);
                sumNon += d;
            }
        }

        reduce(nNon, sumOp<label>());
        reduce(maxNon, maxOp<scalar>());
        reduce(sumNon, sumOp<scalar>());
        sumNon /= returnReduce(D.size(), sumOp<label>());

        InfoIn("faMatrix<Type>::relax(const scalar alpha)")
            << "Matrix dominance test for " << psi_.name() << nl
            << "    number of non-dominant cells   : " << nNon << nl
            << "    maximum relative non-dominance : " << maxNon << nl
            << "    average relative non-dominance : " << sumNon << nl
            << endl;
    }

    // Ensure the matrix is diagonally dominant...
    // Assumes that the central coefficient is positive and ensures it is
    forAll(D, facei)
    {
        D[facei] = max(mag(D[facei]), sumOff[facei]);
    }


    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= component(iCoeffs[face], 0);
                }
            }
            else
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= cmptMin(iCoeffs[face]);
                }
            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.primitiveField();
}


template<class Type>
void faMatrix<Type>::relax()
{
    word name = psi_.select
    (
        psi_.mesh().data::template lookupOrDefault<bool>
        ("finalIteration", false)
    );

    if (psi_.mesh().solution().relaxEquation(name))
    {
        relax(psi_.mesh().solution().equationRelaxationFactor(name));
    }
}

template<class Type>
void Foam::faMatrix<Type>::relax(word name)
{
    if (psi_.mesh().solution().relaxEquation(name))
    {
        relax(psi_.mesh().solution().equationRelaxationFactor(name));
    }
}

template<class Type>
void Foam::faMatrix<Type>::boundaryManipulate
(
    typename GeometricField<Type, faPatchField, areaMesh>::
        Boundary& bFields
)
{
    forAll(bFields, patchI)
    {
        bFields[patchI].manipulateMatrix(*this);
    }
}

template<class Type>
tmp<scalarField> faMatrix<Type>::D() const
{
    tmp<scalarField> tdiag(new scalarField(diag()));
    addCmptAvBoundaryDiag(tdiag.ref());
    return tdiag;
}

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::faMatrix<Type>::DD() const
{
    tmp<Field<Type>> tdiag(pTraits<Type>::one*diag());

    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (!ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchAddr(patchI),
                internalCoeffs_[patchI],
                tdiag()
            );
        }
    }

    return tdiag;
}

template<class Type>
tmp<areaScalarField> faMatrix<Type>::A() const
{
    tmp<areaScalarField> tAphi
    (
        new areaScalarField
        (
            IOobject
            (
                "A("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimArea,
            zeroGradientFaPatchScalarField::typeName,
            psi().boundaryField().patchTypes()
        )
    );

    tAphi.ref().ref() = D()/psi_.mesh().S();
    tAphi.ref().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>> faMatrix<Type>::H() const
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tHphi
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions_/dimArea,
            zeroGradientFaPatchScalarField::typeName,
            psi().boundaryField().patchTypes()
        )
    );

    GeometricField<Type, faPatchField, areaMesh>& Hphi = tHphi.ref();
    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi_.primitiveField().component(cmpt);

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        Hphi.primitiveRefField().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    Hphi.primitiveFieldRef() += lduMatrix::H(psi_.primitiveField()) + source_;
    addBoundarySource(Hphi.primitiveFieldRef(), psi_);

    Hphi.primitiveFieldRef() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}

template<class Type>
Foam::tmp<Foam::areaScalarField> Foam::faMatrix<Type>::H1() const
{
    tmp<areaScalarField> tH1
    (
        new areaScalarField
        (
            IOobject
            (
                "H(1)",
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/(dimArea*psi_.dimensions()),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& H1_ = tH1.ref();

    H1_.primitiveFieldRef() = lduMatrix::H1();

    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                lduAddr().patchAddr(patchI),
                boundaryCoeffs_[patchI].component(0),
                H1_
            );
        }
    }

    H1_.primitiveFieldRef() /= psi_.mesh().S();
    H1_.correctBoundaryConditions();

    return tH1;
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>> faMatrix<Type>::
flux() const
{
    if (!psi_.mesh().schemes().fluxRequired(psi_.name()))
    {
        FatalErrorInFunction
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary of faSchemes"
            << abort(FatalError);
    }

    // construct GeometricField<Type, faePatchField, faEdgeMesh>
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tfieldFlux
    (
        new GeometricField<Type, faePatchField, faEdgeMesh>
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions()
        )
    );
    GeometricField<Type, faePatchField, faEdgeMesh>& fieldFlux = tfieldFlux.ref();

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        fieldFlux.primitiveFieldRef().replace
        (
            cmpt,
            lduMatrix::faceH(psi_.primitiveField().component(cmpt))
        );
    }

    FieldField<Field, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchI)
    {
        InternalContrib[patchI] =
            cmptMultiply
            (
                InternalContrib[patchI],
                psi_.boundaryField()[patchI].patchInternalField()
            );
    }

    FieldField<Field, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchI)
    {
        if (psi_.boundaryField()[patchI].coupled())
        {
            NeighbourContrib[patchI] =
                cmptMultiply
                (
                    NeighbourContrib[patchI],
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
        }
    }

    forAll(fieldFlux.boundaryField(), patchI)
    {
        fieldFlux.boundaryFieldRef()[patchI] =
            InternalContrib[patchI] - NeighbourContrib[patchI];
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void faMatrix<Type>::operator=(const faMatrix<Type>& famv)
{
    if (this == &famv)
    {
        FatalErrorInFunction
            << "attempted to assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(famv.psi_))
    {
        FatalErrorInFunction
            << "different fields"
            << abort(FatalError);
    }

    dimensions_ = famv.dimensions_;
    lduMatrix::operator=(famv);
    source_ = famv.source_;
    internalCoeffs_ = famv.internalCoeffs_;
    boundaryCoeffs_ = famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, faePatchField, faEdgeMesh>
        (*famv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void faMatrix<Type>::operator=(const tmp<faMatrix<Type>>& tfamv)
{
    operator=(tfamv());
    tfamv.clear();
}


template<class Type>
void faMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void faMatrix<Type>::operator+=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ += famv.dimensions_;
    lduMatrix::operator+=(famv);
    source_ += famv.source_;
    internalCoeffs_ += famv.internalCoeffs_;
    boundaryCoeffs_ += famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, faePatchField, faEdgeMesh>
        (
            *famv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void faMatrix<Type>::operator+=(const tmp<faMatrix<Type>>& tfamv)
{
    operator+=(tfamv());
    tfamv.clear();
}


template<class Type>
void faMatrix<Type>::operator-=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ -= famv.dimensions_;
    lduMatrix::operator-=(famv);
    source_ -= famv.source_;
    internalCoeffs_ -= famv.internalCoeffs_;
    boundaryCoeffs_ -= famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, faePatchField, faEdgeMesh>
        (-*famv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void faMatrix<Type>::operator-=(const tmp<faMatrix<Type>>& tfamv)
{
    operator-=(tfamv());
    tfamv.clear();
}


template<class Type>
void faMatrix<Type>::operator+=
(
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().S()*su.primitiveField();
}


template<class Type>
void faMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void faMatrix<Type>::operator-=
(
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().S()*su.primitiveField();
}


template<class Type>
void faMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void faMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= su.mesh().S()*su;
}


template<class Type>
void faMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += su.mesh().S()*su;
}


template<class Type>
void faMatrix<Type>::operator*=
(
    const areaScalarField& vsf
)
{
    dimensions_ *= vsf.dimensions();
    lduMatrix::operator*=(vsf.field());
    source_ *= vsf.field();

    forAll(vsf.boundaryField(), patchI)
    {
        const faPatchScalarField& psf = vsf.boundaryField()[patchI];

        if (psf.coupled())
        {
            internalCoeffs_[patchI] *= psf.patchInternalField();
            boundaryCoeffs_[patchI] *= psf.patchNeighbourField();
        }
        else
        {
            internalCoeffs_[patchI] *= psf.patchInternalField();
            boundaryCoeffs_[patchI] *= psf;
        }
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorInFunction
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void faMatrix<Type>::operator*=
(
    const tmp<areaScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void faMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const faMatrix<Type>& fam1,
    const faMatrix<Type>& fam2,
    const char* op
)
{
    if (&fam1.psi() != &fam2.psi())
    {
        FatalErrorInFunction
            << "incompatible fields for operation "
            << endl << "    "
            << "[" << fam1.psi().name() << "] "
            << op
            << " [" << fam2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fam1.dimensions() != fam2.dimensions())
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam1.psi().name() << fam1.dimensions()/dimArea << " ] "
            << op
            << " [" << fam2.psi().name() << fam2.dimensions()/dimArea << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const faMatrix<Type>& fam,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const char* op
)
{
    if (dimensionSet::debug && fam.dimensions()/dimArea != vf.dimensions())
    {
        FatalErrorInFunction
            <<  "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam.psi().name() << fam.dimensions()/dimArea << " ] "
            << op
            << " [" << vf.name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const faMatrix<Type>& fam,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && fam.dimensions()/dimArea != dt.dimensions())
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam.psi().name() << fam.dimensions()/dimArea << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
solverPerformance solve
(
    faMatrix<Type>& fam,
    Istream& solverControls
)
{
    return fam.solve(solverControls);
}

template<class Type>
solverPerformance solve
(
    const tmp<faMatrix<Type>>& tfam,
    Istream& solverControls
)
{
    solverPerformance solverPerf =
        const_cast<faMatrix<Type>&>(tfam()).solve(solverControls);

    tfam.clear();
    return solverPerf;
}


template<class Type>
solverPerformance solve(faMatrix<Type>& fam)
{
    return fam.solve();
}

template<class Type>
solverPerformance solve(const tmp<faMatrix<Type>>& tfam)
{
    solverPerformance solverPerf =
        const_cast<faMatrix<Type>&>(tfam()).solve();

    tfam.clear();
    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>> operator+
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() += B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<faMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() -= B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<faMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator==
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
tmp<faMatrix<Type>> operator==
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
tmp<faMatrix<Type>> operator==
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
tmp<faMatrix<Type>> operator==
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator+
(
    const tmp<faMatrix<Type>>& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator+
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator+
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator+
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<faMatrix<Type>>& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator-
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.mesh().S()*su.primitiveField();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC.ref().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator+
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC.ref().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += su.value()*tC.ref().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.value()*tC.ref().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator-
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.value()*tC.ref().psi().mesh().S();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator==
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator==
(
    const tmp<faMatrix<Type>>& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().S()*su.primitiveField();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator==
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator==
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().primitiveField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator==
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += A.psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator==
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tC.ref().psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator*
(
    const areaScalarField& vsf,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() *= vsf;
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator*
(
    const tmp<areaScalarField>& tvsf,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() *= tvsf;
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator*
(
    const areaScalarField& vsf,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= vsf;
    return tC;
}

template<class Type>
tmp<faMatrix<Type>> operator*
(
    const tmp<areaScalarField>& tvsf,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tvsf;
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator*
(
    const dimensioned<scalar>& ds,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() *= ds;
    return tC;
}


template<class Type>
tmp<faMatrix<Type>> operator*
(
    const dimensioned<scalar>& ds,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= ds;
    return tC;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const faMatrix<Type>& fam)
{
    os  << static_cast<const lduMatrix&>(fam) << nl
        << fam.dimensions_ << nl
        << fam.source_ << nl
        << fam.internalCoeffs_ << nl
        << fam.boundaryCoeffs_ << endl;

    os.check("Ostream& operator<<(Ostream&, faMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "faMatrices/faMatrix/faMatrixSolve.C"

// ************************************************************************* //
