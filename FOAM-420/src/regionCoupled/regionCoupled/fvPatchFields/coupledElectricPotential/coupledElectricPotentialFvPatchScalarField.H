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
    (c) 2010-2019 Esi Ltd.
    (c) 2012-2013 OpenFOAM Foundation

Class
    Foam::coupledElectricPotentialFvPatchScalarField

Description
    Region coupled implicit boundary condition for electric potential
    supporting optional contact resistance. The fvPatch is treated as uncoupled
    from the delta point of view, as all regions are not always available for
    interpolations etc. In the mesh the fvPatch is an interface and is
    incorporated into the matrix implicitly.
    The electrical conductivity, \c sigma, is retrieved from the
    mesh database using a lookup.
    \heading Patch usage
    \table
        Property     | Description             | Required    | Default value
        contactResistance | Total electrical contact resistance [Ohm] | no |
        sigmaLayers  | List of electrical conductivites per layer [S/m] | no |
        thicknessLayers | List of thicknesses per layer [m] | no |
    \endtable
    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            coupledElectricPotential;
        sigmaLayers     ( 1e3 1e5 );
        thicknessLayers (0.001 0.01 );
        value           uniform 0;
    }
    \endverbatim
    Needs to be on underlying mappedWallFvPatch.

SourceFiles
    coupledElectricPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledElectricPotentialFvPatchScalarField_H
#define coupledElectricPotentialFvPatchScalarField_H

#include "fvMesh/fvPatches/derived/mapped/mappedWallFvPatch.H"
#include "matrices/LduMatrix/LduMatrix/LduInterfaceField/LduInterfaceField.H"
#include "fvPatchFields/regionCoupled/regionCoupledFvPatchField.H"
#include "fvPatchFields/electricalBoundaryBase/electricalBoundaryBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class coupledElectricPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class coupledElectricPotentialFvPatchScalarField
:
    public regionCoupledFvPatchField<scalar>,
    public electricalBoundaryBase
{
private:

    // Private data

    //- Owner cell weight
    scalarField lWeights_;

    //- Explicit correction to face value
    scalarField fCorr_;

    //- Total contact resistance
    scalar contactResistance_;

    //- Conductivity of layers
    scalarList sigmaLayers_;

    //- Thickness of layers
    scalarList thicknessLayers_;

    //- Inverse of contact conductance per unit area (computed from above)
    scalar contactInsulance_;


    // Private functions

    //- Local weight for this coupled field
    tmp<scalarField> weights() const;


public:

    //- Runtime type information
    TypeName("coupledElectricPotential");


    // Constructors

    //- Construct from patch and internal field
    coupledElectricPotentialFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    coupledElectricPotentialFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given coupledElectricPotentialFvPatchScalarField
    // onto a new patch
    coupledElectricPotentialFvPatchScalarField
    (
        const coupledElectricPotentialFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    coupledElectricPotentialFvPatchScalarField
    (
        const coupledElectricPotentialFvPatchScalarField&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<scalar>> clone() const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new coupledElectricPotentialFvPatchScalarField(*this)
            );
    }

    //- Construct as copy setting internal field reference
    coupledElectricPotentialFvPatchScalarField
    (
        const coupledElectricPotentialFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const override
    {
        return tmp<fvPatchField<scalar>>
            (
                new coupledElectricPotentialFvPatchScalarField(*this, iF)
            );
    }


    //- Destructor
    virtual ~coupledElectricPotentialFvPatchScalarField()
    {}


    // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        ) override;

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchScalarField&,
            const labelList&
        ) override;


    // Evaluation functions

        //- Return patch-normal gradient
        virtual tmp<scalarField> snGrad() const override;

        //- Return patch-normal gradient
        //  Note: the deltaCoeffs supplied are not used
        virtual tmp<scalarField> snGrad
        (
            const scalarField& deltaCoeffs
        ) const override;


        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        ) override;

        virtual void updateCoeffs() override;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField
        virtual tmp<Field<scalar>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField
        virtual tmp<Field<scalar>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const override;

        //- Return the explicit correction region boundary face value
        virtual tmp<Field<scalar>> valueCorrection() const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientInternalCoeffs() const override;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const override;

        //- Return the explicit correction of the normal gradient
        virtual tmp<Field<scalar>> gradientCorrection() const;

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
            return tmp<scalarField>(fCorr_);
        }

        //- Return explicit correction to current boundary face value
        tmp<scalarField> faceCorrEval() const override
        {
            return tmp<scalarField>(fCorr_);
        }


    //- Write
    virtual void write(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
