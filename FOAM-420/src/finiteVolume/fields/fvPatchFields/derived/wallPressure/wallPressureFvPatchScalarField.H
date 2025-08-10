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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::wallPressureFvPatchScalarField

Group

Description
    Pressure boundary for single phase incompressible isothermal flows.
    Supplies extrapolated wall pressure gradient to the momentum equation,
    making the source term second order.
    Uses zero gradient for implicit Laplacian terms to ensure zero wall flux.


SeeAlso
    Foam::fixedGradientFvPatchField

SourceFiles
    wallPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallPressureFvPatchScalarFields_H
#define wallPressureFvPatchScalarFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchFields.H"
#include "primitives/bools/Switch/Switch.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class wallPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallPressureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Current local index (used for updating)
        label locInterpIndex_;

        //- addressing for cells with multiple faces on non-wallPressure patches
        autoPtr<List<List<Tuple2<label, label>>>> multiFaceAddrPtr_;

        //- addressing for cells with faces on wallPressure patches
        autoPtr<List<List<Tuple2<label, label>>>> wallPressureFaceAddrPtr_;

        //- start angle for damping gradient correction in skew boundary cells
        scalar startDampingAngle_;

        //- non-orthogonality angle at which the  gradient becomes zero
        scalar zeroGradientAngle_;

        //- gradient reduction factor (for reuse)
        scalarField alphaG_;

        //- gradient relaxation factor
        scalar relaxG_;

        //- update addressing
        bool updateAddressing_;

        //- Global update index for interpolated face field
        static label globInterpIndex;

        //- Interpolated face field
        static autoPtr<surfaceScalarField> facePressure;

    // Private member functions

        //- calculate addressing for cells on boundary with multiple faces
        // on wallPressure patches
        void calculateMultiFaceAdressing();

        //- calculate face normal cell gradient implicitly
        scalar faceNormalGradient
        (
            const fvPatchScalarField& cPatch,
            label faceI,
            const vector& nf,
            const List<Tuple2<label, label>>& otherCellBoundaryFaces,
            label cellI
        ) ;

        //- calculate explicit cell gradient using Gauss
        vector cellGradient
        (
            label cellI,
            const Tuple2<label, label>& startCellBF,
            const List<Tuple2<label, label>>& wallPressureBFs,
            const List<Tuple2<label, label>>& otherCellBFs
        ) const;

        //- implicit update boundary gradient via extrapolation
        void extrapolateGradient();

public:

    //- Runtime type information
    TypeName("wallPressure");


    // Constructors

        //- Construct from patch and internal field
        wallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wallPressureFvPatchScalarField onto
        //  a new patch
        wallPressureFvPatchScalarField
        (
            const wallPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallPressureFvPatchScalarField
        (
            const wallPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new wallPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wallPressureFvPatchScalarField
        (
            const wallPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new wallPressureFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    ~wallPressureFvPatchScalarField();

    // Member functions

        //- return gradient reduction factor
        inline const scalarField& gradReductionFactor() const
        {
            return alphaG_;
        }

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<scalarField> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<scalarField> gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
