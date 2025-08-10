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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2010-2022 Esi Ltd.

Class
    Foam::regionCoupledFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    Abstract base class for region-coupled patch fields.

SourceFiles
    regionCoupledFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledFvPatchField_H
#define regionCoupledFvPatchField_H

#include "fields/CoeffField/CoeffField.H"
#include "fields/fvPatchFields/basic/coupled/coupledFvPatchField.H"
#include "fvMesh/fvPatches/derived/mapped/mappedFvPatch.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coupledFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class regionCoupledFvPatchField
:
    public coupledFvPatchField<Type>
{
    //- Check for if the boundary is next to inactive region
    bool lookupInactiveRegions();

    //- Check for if the boundary is next to inactive region
    fvPatchField<Type>* setHeFromT
    (
        const fvPatch& p,
        const DimensionedField<Type, volMesh>& iF
    );

protected:

    // Protected data

        //- Local reference to region-coupled patch
        const mappedFvPatch& regionCoupledPatch_;

        //- Is it region coupled boundary
        // (not if the other region ins't being solved)
        bool isRegionCoupled_;

        //- Name of coupled field - defaults to this field name
        mutable word neighbourFieldName_;

        //- Name of coupled group name - defaults to this one
        word neighbourPhaseName_;

        //- Name of coupled region - if unspecified, automatically selects
        //  the region with a corresponding region-coupled patch field, if
        //  there is only one
        word neighbourRegion_;

        //- Optional boundary for switch-off region
        autoPtr<fvPatchField<Type>> nonCoupledRegionBoundary_;

public:

    //- Runtime type information
    TypeName("regionCoupled");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        regionCoupledFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given regionCoupledFvPatchField onto a new patch
        regionCoupledFvPatchField
        (
            const regionCoupledFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        regionCoupledFvPatchField
        (
            const regionCoupledFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        regionCoupledFvPatchField
        (
            const regionCoupledFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const override
        {
            return tmp<fvPatchField<Type>>
            (
                new regionCoupledFvPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchField<Type>>
            (
                new regionCoupledFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Do not attempt interpolations across region couplings
            virtual bool coupled() const override
            {
                return false;
            }

            virtual bool regionCoupled() const override
            {
                return isRegionCoupled_;
            }

            //- Return name of other side field that is coupled
            virtual word coupledField() const;

            //- Return neighbour field of internal field
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            ) const override;

            //- Return object registry in other region
            const objectRegistry& neighbourDb() const;

            //- Return coupled patch field in other region
            virtual const regionCoupledFvPatchField<Type>&
            neighbourFvPatchField() const;

            const word& sampleRegion() const
            {
                return regionCoupledPatch_.sampleRegion();
            }

            //- Inactive boundary
            autoPtr<fvPatchField<Type>>& nonCoupledRegionBoundary()
            {
                return nonCoupledRegionBoundary_;
            }

            //- Inactive boundary
            const autoPtr<fvPatchField<Type>>&
            nonCoupledRegionBoundary() const
            {
                return nonCoupledRegionBoundary_;
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            ) override;

            //- Coupled version should not be used
            virtual tmp<Field<Type>> snGrad
            (
                const scalarField& deltaCoeffs
            ) const override;

            //- Coupled version should not be used
            virtual tmp<Field<Type>> snGrad() const override;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const override;

            //- Coupled version should not be used
            virtual tmp<Field<Type>> gradientInternalCoeffs
            (
                const scalarField& deltaCoeffs
            ) const override;

            //- Return the matrix coupled coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const override;

            //- Coupled version should not be used
            virtual tmp<Field<Type>> gradientBoundaryCoeffs
            (
                const scalarField& deltaCoeffs
            ) const override;

            //- Return patch value for explicit gradient evaluation
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const override;


        // Coupled interface functionality

            //- Convert boundary/internal coeffs to region-coupled
            //  boundary/internal coeffs (from coeffs of boundary face value and
            //  internal value to coeffs of neighbour & internal value)
            virtual void regionCoupledBoundaryCoeffs
            (
                const fvMatrix<Type>& matrix,
                const Field<Type>& bouCoeffs,
                const Field<Type>& intCoeffs,
                Field<Type>& coupledBouCoeffs,
                Field<Type>& coupledIntCoeffs
            )
            {}

            //- Return new explicit correction to boundary face value
            virtual tmp<Field<Type>> faceCorr() const
            {
                return
                    tmp<Field<Type>>
                    (
                        new Field<Type>(this->size(), pTraits<Type>::zero)
                    );
            }

            //- Return explicit correction to current boundary face value
            virtual tmp<Field<Type>> faceCorrEval() const
            {
                return
                    tmp<Field<Type>>
                    (
                        new Field<Type>(this->size(), pTraits<Type>::zero)
                    );
            }

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction,
                const Pstream::commsTypes commsType
            ) const override;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const override;

            //- Return the interface type
            virtual const word& interfaceFieldType() const override
            {
                return regionCoupledPatch_.regionCoupleType();
            }

            //- Interpolate from neighbour using region coupling patch
            virtual tmp<scalarField> interpolateFromNeighbour
            (
                const scalarField&,
                const UList<scalar>& defaultValues,
                const direction cmpt
            ) const;

            //- Interpolate from neighbour using region coupling patch
            virtual tmp<Field<Type>> interpolateFromNeighbour
            (
                const Field<Type>&, const UList<Type>& defaultValues
            ) const;

            //- Interpolate to neighbour using region coupling patch
            virtual tmp<Field<Type>> interpolateToNeighbour
            (
                const Field<Type>&, const UList<Type>& defaultValues
            ) const;

        // I/O

            virtual void write(Ostream& os) const override;
};

// Specialisations
template<>
fvPatchField<scalar>* regionCoupledFvPatchField<scalar>::setHeFromT
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
