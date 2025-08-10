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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/schemes/limitedLinearUpwind/limitedLinearUpwind.H"
#include "finiteVolume/fvc/fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceInterpolationTypeScheme(limitedLinearUpwind, scalar);
    makeSurfaceInterpolationTypeScheme(limitedLinearUpwind, vector);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::limitedLinearUpwind<Type>::fLimiter
(
    const scalar& xRaw
) const
{
    scalar x = xRaw;

    //-bound x
    x = min(x, 2);
    x = max(x, 0);

    scalar fx = 1.0;
    if (x==0)
    {
        fx = 0.0;
    }
    else if (x==2)
    {
        fx = 1.0;
    }
    else
    {
        if (limiterTypeName_ == ("step"))
        {
            //-original BJ
            fx = min(x, 1.0);
        }
        else if
        (
            limiterTypeName_ == ("functional")
         || limiterTypeName_ == ("default")
        )
        {
            //- Venkatakrishnana
            fx = (x*x + 2.0*x)/(x*x + x + 2.0);
        }
        else if (limiterTypeName_ == ("R3"))
        {
            //- Check [8]
            fx = (pow(x,3) + 4*x)/(pow(x,3) + x*x + 4);
        }
        else if (limiterTypeName_ == ("R4"))
        {
            //- Check [8]
            fx = (pow(x,4) + 2*pow(x,3) - 4*pow(x,2) + 8*x)
               /(pow(x,4) + pow(x,3) + 2*pow(x,2) - 4*x + 8);
        }
        else if (limiterTypeName_ == ("R5"))
        {
            //- Check [8]
            fx = (pow(x,5) + 8* pow(x,3) - 16*pow(x,2) + 16*x)
               /((pow(x,5) + pow(x,4) + 8*pow(x,2) - 16*x +16));
        }
    }

    //-sanity bounding. This should be ok from the above
    fx = min(fx, 1.0);
    fx = max(fx, 0.0);

    return fx;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::limitedLinearUpwind<Type>::computePsiLimiter
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
) const
{

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tpsi
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "psiLimiter_"+vf.name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh(),
            dimensioned<Type>("psi", dimless, pTraits<Type>::one)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& limiterField = tpsi.ref();
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr();

    //- Compute min-max value of nei cells in each direction
    tmp<Field<Type>> tminLocalCellNeiValue(cellNeiOp(vf, minOp<scalar>()));
    tmp<Field<Type>> tmaxLocalCellNeiValue(cellNeiOp(vf, maxOp<scalar>()));
    const Field<Type>& minLocalCellNeiValue = tminLocalCellNeiValue();
    const Field<Type>& maxLocalCellNeiValue = tmaxLocalCellNeiValue();

    //- Compute Gradient
    tmp<fv::gradScheme<Type>> gradScheme
    (
        fv::gradScheme<Type>::New
        (
            this->mesh(),
            vf.db(),
            this->mesh().schemes().gradScheme("grad(" + vf.name() + ')')
        )
    );

    typedef
    GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        fvPatchField,
        volMesh
    > GradFieldType;

    tmp<GradFieldType> tgradV = gradScheme().grad(vf);

    const GradFieldType& gradVf = tgradV();

    //- Calculate denominator --> \nabla v_P d and limiter
    const volVectorField& C = this->mesh().C();
    const surfaceVectorField& Cf = this->mesh().Cf();
    const Field<Type>& sfCorrInt = sfCorr.primitiveField();

    const vectorField& CfInt = Cf.internalField();
    const vectorField& CInt = C.internalField();

    forAll(sfCorrInt, fI)
    {
        const labelList& owner = this->mesh().owner();
        const labelList& neighbour = this->mesh().neighbour();

        // Investigate face-based limiter in both sides to take all possible
        // situations
        //- owner
        {
            const label cI = owner[fI];
            Type denominatorI = (CfInt[fI] - CInt[cI])&gradVf[cI];

            psiLimiterType
            (
                denominatorI,
                vf[cI],
                maxLocalCellNeiValue[cI],
                minLocalCellNeiValue[cI],
                limiterField[fI]
            );
        }
        //- neighbour
        {
            const label cI = neighbour[fI];
            Type denominatorI = (CfInt[fI] - CInt[cI])&gradVf[cI];

            psiLimiterType
            (
                denominatorI,
                vf[cI],
                maxLocalCellNeiValue[cI],
                minLocalCellNeiValue[cI],
                limiterField[fI]
            );
        }
    }
    //- coupled boundaries
    forAll(sfCorr.boundaryField(), pI)
    {
        const fvPatchField<Type>& vfP = vf.boundaryField()[pI];

        if (vfP.coupled())
        {
            const labelList& pFaceCells = vfP.patch().faceCells();
            const Field<Type> denominatorI = sfCorr.boundaryField()[pI];
            const vectorField& pCf = Cf.boundaryField()[pI];
            Field<Type>& bLimiterField = limiterField.boundaryFieldRef()[pI];

            forAll(pFaceCells, fI)
            {
                const label cI = pFaceCells[fI];
                Type denominatorI = (pCf[fI] - CInt[cI])&gradVf[cI];
                psiLimiterType
                (
                    denominatorI,
                    vf[cI],
                    maxLocalCellNeiValue[cI],
                    minLocalCellNeiValue[cI],
                    bLimiterField[fI]
                );
            }
        }
    }

    return tpsi;
}


template<class Type>
template<class BinaryOp>
Foam::tmp<Foam::Field<Type>> Foam::limitedLinearUpwind<Type>::cellNeiOp
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const BinaryOp& bop
) const
{
    //-Loop in cells and do the bop operation. Clean way to take the min-max in
    //cell-cells. Works in a component-by-component basis

    tmp<Field<Type>> tField
    (
        new Field<Type>(vf.primitiveField())
    );
    Field<Type>& field = tField.ref();

    const labelList& owner = this->mesh().owner();
    const labelList& neighbour = this->mesh().neighbour();
    const Field<Type>& vfInt = vf.internalField();

    forAll(neighbour, fI)
    {
        label own = owner[fI];
        label nei = neighbour[fI];

        Type& fOwn = field[own];
        Type& fNei = field[nei];

        const Type& vNei = vfInt[nei];
        const Type& vOwn = vfInt[own];

        bopType(fOwn, vNei, bop);
        bopType(fNei, vOwn, bop);
    }

    forAll(vf.boundaryField(), pI)
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[pI];
        const labelList& pFaceCells = pvf.patch().faceCells();

        if (pvf.coupled())
        {
            const Field<Type> pvfNei(pvf.patchNeighbourField());

            forAll(pFaceCells, fI)
            {
                const label cI = pFaceCells[fI];
                Type& fOwn = field[cI];
                const Type& vNei = pvfNei[fI];
                bopType(fOwn, vNei, bop);
            }
        }
        else
        {
            forAll(pFaceCells, fI)
            {
                label cI = pFaceCells[fI];
                Type& fOwn = field[cI];
                const Type& vfOwn = pvf[fI];
                bopType(fOwn, vfOwn, bop);
            }
        }
    }

    return tField;
}

template<class Type>
void Foam::limitedLinearUpwind<Type>::psiLimiter
(
    const scalar& denominatorI,
    const scalar& vfI,
    const scalar& maxI,
    const scalar& minI,
    scalar& psiI
) const
{
    //- Init with no limiter
    scalar x = 2.0;

    //- Check both sides min/max based on flux sign
    if (denominatorI>0)
    {
        //- division by very small number. Small addition
        scalar denpSmall = denominatorI+SMALL;

        x = (maxI - vfI)/denpSmall;
        x = min(x, -(minI - vfI)/denpSmall);
    }
    else if (denominatorI<0)
    {
        //- division by very small number. Small addition
        scalar denmSmall = denominatorI-SMALL;

        x = -(maxI - vfI)/denmSmall;
        x = min(x, (minI - vfI)/denmSmall);
    }
    else
    {
        //-leave x = 2 --> no limiter
    }

    //- calculate psi = f(x)
    //- min because a face can be investigated twice (own-nei)
    psiI = min(fLimiter(x), psiI);
}


template<class Type>
void Foam::limitedLinearUpwind<Type>::psiLimiterType
(
    const Type& denominatorI,
    const Type& vfI,
    const Type& maxI,
    const Type& minI,
    Type& psiI
) const
{
    const direction nCmpts = pTraits<Type>::nComponents;
    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        const scalar& denominatorICmpt = component(denominatorI, cmptI);

        const scalar& vfCmpt = component(vfI, cmptI);
        const scalar& maxCmpt = component(maxI, cmptI);
        const scalar& minCmpt = component(minI, cmptI);
        scalar& limitCmpt = setComponent(psiI, cmptI);

        psiLimiter
        (
            denominatorICmpt, vfCmpt, maxCmpt, minCmpt, limitCmpt
        );
    }
}


template<class Type>
void Foam::limitedLinearUpwind<Type>::writeFieldsforDebug
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tpsi
) const
{
    if (this->mesh().time().outputTime())
    {
        tpsi->write();

        GeometricField<Type, fvPatchField, volMesh> maxCellLim
        (
            IOobject
            (
                "cell"+tpsi->name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh(),
            dimensioned<Type>("psiCell", dimless, pTraits<Type>::zero)
        );
        const labelList& owner = this->mesh().owner();
        const labelList& neighbour = this->mesh().neighbour();
        surfaceScalarField magSf = this->mesh().magSf();
        scalarField sumSF (this->mesh().cells().size(), 0);
        forAll(neighbour, fI)
        {
            label own = owner[fI];
            label nei = neighbour[fI];
            maxCellLim.ref()[own] +=
                tpsi().primitiveField()[fI]*magSf[fI];
            maxCellLim.ref()[nei] +=
                tpsi().primitiveField()[fI]*magSf[fI];

            sumSF[own] += magSf[fI];
            sumSF[nei] += magSf[fI];
        }
        forAll(maxCellLim, cI)
        {
            maxCellLim.ref()[cI] /= sumSF[cI];
        }
        maxCellLim.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::limitedLinearUpwind<Type>::limitedLinearUpwind
(
    const fvMesh& mesh,
    const objectRegistry& db,
    Istream& schemeData
)
:
    linearUpwind<Type>(mesh, db, schemeData),
    limiterTypeName_(schemeData)
{}


template<class Type>
Foam::limitedLinearUpwind<Type>::limitedLinearUpwind
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
:
    linearUpwind<Type>(mesh, faceFlux, schemeData),
    limiterTypeName_(schemeData)
{}

// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::limitedLinearUpwind<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        tlinearUpwindCorrection = linearUpwind<Type>::correction(vf);

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tpsi =
        computePsiLimiter(vf, tlinearUpwindCorrection);

    if (debug)
    {
        writeFieldsforDebug(tpsi);
    }

    return cmptMultiply(tpsi(), tlinearUpwindCorrection());
}

// ************************************************************************* //
