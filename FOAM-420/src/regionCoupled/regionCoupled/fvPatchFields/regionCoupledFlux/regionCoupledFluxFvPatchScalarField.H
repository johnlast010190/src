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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2010-2021 Esi Ltd.

Class
    Foam::regionCoupledFluxFvPatchScalarField

Description
    Region coupled implicit boundary condition based on equation flux equality.
    Supports boundary sources, and optional thin contact resistance.
    The fvPatch is treated as uncoupled from the delta point of view, as
    all regions are not always available for interpolations etc.
    In the mesh the fvPatch is an interface and is incorporated
    into the matrix implicitly.
    \heading Patch usage
    \table
        Property     | Description             | Required    | Default value
        contactResistance | Total contact resistance | no |
        conductivityLayers  | List of conductivity per layer | no |
        thicknessLayers | List of thickness per layer | no |
    \endtable
    Only one of either contactResistance or conductivityLayers and
    thicknesslayers may be specified.
    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            regionCoupledFlux;
        conductivityLayers   ( 1e3 1e5 );
        thicknessLayers     (0.001 0.01 );
        value           uniform 300;
    }
    \endverbatim
    Needs to be on underlying mappedFvPatch.

SourceFiles
    regionCoupledFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledFluxFvPatchScalarField_H
#define regionCoupledFluxFvPatchScalarField_H

#include "fvMesh/fvPatches/derived/mapped/mappedWallFvPatch.H"
#include "matrices/LduMatrix/LduMatrix/LduInterfaceField/LduInterfaceField.H"
#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class regionCoupledFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledFluxFvPatchScalarField
:
    public regionCoupledFvPatchField<scalar>
{
protected:

    // Protected data

        //- Stored intermediates
        scalarField Cb_;
        scalarField Ci_;
        scalarField fluxCorr_;

        //- Stored for reconstruction of face value
        scalarField faceBCoeff_;
        scalarField faceICoeff_;
        scalarField faceCorr_;

        //- Total contact resistance
        scalar contactResistance_;

        //- Diffusivity of layers
        scalarList conductivityLayers_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Inverse of contact conductivity per unit area (computed from either
        //  layers or contact resistance above)
        scalar contactInsulance_;

        //- Net patch flux source (used for damping)
        //  The value stored in QsumNbr_ is actually used in the neighbour's
        //  calculation and vice-versa, so that the field/patch sizing matches
        scalarField Qsum_;
        scalarField QsumNbr_;

        //- Heat source under-relaxation factor
        scalar qURF_;


    // Protected functions

        //- Overridable functions to scale and offset such that
        //  fieldValue = scaleFactor*(matchedValue+offset)
        //  where matchedValue is continuous across the interface
        virtual tmp<scalarField> fieldScale(const tmp<scalarField>& f) const
        {
            return tmp<scalarField>(f.ptr());
        }

        virtual tmp<scalarField> invFieldScale(const tmp<scalarField>& f) const
        {
            return tmp<scalarField>(f.ptr());
        }

        virtual tmp<scalarField> fieldOffset(const tmp<scalarField>& f) const
        {
            return tmp<scalarField>(f.ptr());
        }

        virtual tmp<scalarField> negFieldOffset(const tmp<scalarField>& f) const
        {
            return tmp<scalarField>(f.ptr());
        }

        //- Overridable function to scale the flux such that fluxScale*flux
        //  is equal on either side of the interface
        virtual tmp<scalarField> fluxScale(const tmp<scalarField>& f) const
        {
            return tmp<scalarField>(f.ptr());
        }

        virtual tmp<scalarField> invFluxScale(const tmp<scalarField>& f) const
        {
            return tmp<scalarField>(f.ptr());
        }

public:

    //- Runtime type information
    TypeName("regionCoupledFlux");


    // Constructors

    //- Construct from patch and internal field
    regionCoupledFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    regionCoupledFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given regionCoupledFluxFvPatchScalarField
    // onto a new patch
    regionCoupledFluxFvPatchScalarField
    (
        const regionCoupledFluxFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    regionCoupledFluxFvPatchScalarField
    (
        const regionCoupledFluxFvPatchScalarField&
    );

    //- Construct and return a clone
    tmp<fvPatchField<scalar>> clone() const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new regionCoupledFluxFvPatchScalarField(*this)
            );
    }

    //- Construct as copy setting internal field reference
    regionCoupledFluxFvPatchScalarField
    (
        const regionCoupledFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new regionCoupledFluxFvPatchScalarField(*this, iF)
            );
    }


    //- Destructor
    virtual ~regionCoupledFluxFvPatchScalarField()
    {}


    // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        void autoMap
        (
            const fvPatchFieldMapper&
        ) override;

        //- Reverse map the given fvPatchField onto this fvPatchField
        void rmap
        (
            const fvPatchScalarField&,
            const labelList&
        ) override;

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMapGIB(const gibFvPatchFieldMapper&) override;


    // Evaluation functions

        //- Return patch-normal gradient
        tmp<scalarField> snGrad() const override;

        //- Return patch-normal gradient
        //  Note: the deltaCoeffs supplied are not used
        tmp<scalarField> snGrad
        (
            const scalarField& deltaCoeffs
        ) const override;

        //- Evaluate the patch field
        void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        ) override;

        void updateCoeffs() override;

        //- Return the matrix coefficients of the neighbouring cell value
        //  corresponding to the evaluation of the value of this patchField
        tmp<Field<scalar>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField
        tmp<Field<scalar>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        tmp<Field<scalar>> gradientInternalCoeffs() const override;

        //- Return the matrix coefficients of the neighbouring cell value
        //  corresponding to the evaluation of the gradient of this patchField
        tmp<Field<scalar>> gradientBoundaryCoeffs() const override;

        void manipulateMatrix(fvMatrix<scalar>& matrix) override;

        //- Return the boundary coeffs and explicit correction for use in the
        // region-coupled super-matrix
        void regionCoupledBoundaryCoeffs
        (
            const fvMatrix<scalar>& matrix,
            const scalarField& bouCoeffs,
            const scalarField& intCoeffs,
            scalarField& coupledBouCoeffs,
            scalarField& coupledIntCoeffs
        ) override;

        //- Return explicit correction to new boundary face value
        tmp<scalarField> faceCorr() const override
        {
            return tmp<scalarField>(faceCorr_);
        }

        //- Return explicit correction to current boundary face value
        tmp<scalarField> faceCorrEval() const override
        {
            return tmp<scalarField>(faceCorr_);
        }


    //- Write
    virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
