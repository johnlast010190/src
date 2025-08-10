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

Class
    Foam::limitedLinearUpwind

Group
    grpFvSurfaceInterpolationSchemes

Description
    Scheme based on FROMM scheme which in OpenFOAM is linearUpwind.
    In this scheme, the correction is limited based on
    Barth and Jesperson scheme [1, 2, 3, 4].
    (default) are supported.
    The limiter controls the percentage of the second-order scheme added.
    Thus, the limitedLinearUpwind can be either first, second order or something
    between.

    This scheme follows the Speikreijse monotonicity criterion [2], where the
    values of the reconstructed value should be bounded from the min-max values
    of the neighbor cell elements.

    phi_reconstructed = phi_p + psi \nabla phi_P d

           f((max(phi_N, phi_P) - phi_P) / \nabla phi_P d),  if phi_f > phi_P
    psi =  f((phi_P - min(phi_N, phi_P)) / \nabla phi_P d),  if phi_f < phi_P
           1, if phi_f = phi_P

    f(x) = min(x, 1), original BJ limiter (step)
           (x^2 +2x)/(x^2+x+2), Venkatakrishnan [3] (functional)

    The default is the functional.

    Steps:
    1. The min max of the cell neis are computed
    2. Compute the psi based on the above
    3. Compute the limiter which will have value between 0(upwind) and
       1(linearUpwind).

References

    [1] T. Barth, D.C. Jespersen:
       "The design and application of upwind schemes on unstructured meshes",
       AIAA paper 89-0366, 1989.

    [2] S. Spekreijse:
       "Multigrid solution of monotone, second order discretization of
       hyperbolic conservation laws",
       Mathematics of Computations 45 (1987) 15–21.

    [3] V. Venkatakrishnana:
       "Convergence to steady state solutions of the Euler equations on
       unstructured grids with limiters",
       Journal of Computational Physics 118 (1995) 120–130.

    [4] M. Darwish, F. Moukalled:
       "TVD schemes for unstructured grids",
       International Journal of Heat and Mass Transfer,
       Volume 46, Issue 4, February 2003, Pages 599-611

    [5] P.K. Sweby:
       "High resolution schemes using flux limiters for hyperbolic conservation
       laws." SIAM journal on numerical analysis 21.5, 995-1011, 1984

    [6] B.P. Leonard:
       "Simple high-accuracy resolution program for convective modeling of
       discontinuities", International Journal Numerical Methods
       and Fluids, 8, 1291–1318, 1988.

    [7] H. Gaskell and A.K.C. Lau:
        "Curvature-compensated convective transport: SMART, a new
        boundedness-pre-serving transport algorithm",
        International Journal Numerical Methods and Fluids, 8, 617-641, 1988.

    [8] H. Nishikawa
        "New Unstructured-Grid Limiter Functions",
        National Institute of Aerospace, Hampton, VA 23666, USA

Usage
    div(flux,field)
        bounded Gauss limitedLinearUpwind <gradScheme> <slopeType>;

SourceFiles
    limitedLinearUpwind.C

\*---------------------------------------------------------------------------*/

#ifndef limitedLinearUpwind_H
#define limitedLinearUpwind_H

#include "interpolation/surfaceInterpolation/schemes/linearUpwind/linearUpwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class limitedLinearUpwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedLinearUpwind
:
    public linearUpwind<Type>
{
    // Private Members

        //- limiter name
        const word limiterTypeName_;

    // Private Member Functions

        //- compute f function
        scalar fLimiter(const scalar&) const;

        //- Compute the psi limiter
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        computePsiLimiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        ) const;

        //-
        template<class BinaryOp>
        tmp<Field<Type>> cellNeiOp
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const BinaryOp&
        ) const;

        //- binary operator for type
        template<class BinaryOp>
        void bopType
        (
            Type&,
            const Type&,
            const BinaryOp&
        ) const
        {
            FatalErrorInFunction
              << "Only scalars and vectors are permitted"
              << exit(FatalError);
        }

        //- compute psi
        void psiLimiter
        (
            const scalar&,
            const scalar&,
            const scalar&,
            const scalar&,
            scalar&
        ) const;

        //- compute psi wrapper
        void psiLimiterType
        (
            const Type&,
            const Type&,
            const Type&,
            const Type&,
            Type&
        ) const;

        //- debug helper
        void writeFieldsforDebug
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        ) const;

        //- Disallow default bitwise copy construct
        limitedLinearUpwind(const limitedLinearUpwind&);

        //- Disallow default bitwise assignment
        void operator=(const limitedLinearUpwind&);


public:

    //- Runtime type information
    TypeName("limitedLinearUpwind");


    // Constructors

        //- Construct from mesh and Istream
        limitedLinearUpwind
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        );

        //- Construct from mesh, faceFlux and Istream
        limitedLinearUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return upwind<Type>::weights();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};

// Specialisation for scalars
template<>
template<class BinaryOp>
void Foam::limitedLinearUpwind<scalar>::bopType
(
    scalar& f,
    const scalar& v,
    const BinaryOp& bop
) const
{
    f = bop(f, v);
}


// Specialisation for vectors
template<>
template<class BinaryOp>
void Foam::limitedLinearUpwind<vector>::bopType
(
    vector& f,
    const vector& v,
    const BinaryOp& bop
) const
{
    const direction nCmpts = pTraits<vector>::nComponents;
    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        f.component(cmptI) =
            bop
            (
                f.component(cmptI),
                v.component(cmptI)
            );
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
