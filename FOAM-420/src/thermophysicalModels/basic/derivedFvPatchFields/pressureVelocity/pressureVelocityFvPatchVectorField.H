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
    Foam::pressureVelocityFvPatchVectorField

Description
    pressureVelocity is used in conjunction with the pressureVelocity BC.
    It is used when there is coupling of terms in the expression of the BC,
    like totalPressure and porous BCs.
    In these situations, the pressure is fixed type and the velocity is
    zeroGradient for the normal component and fixedValue for the user-specified
    direction.

Usage

    Example of the boundary condition specification:

    \verbatim
    outlet
    {
        type            pressureVelocity;
        velocityDirection (0.5 0.1 0);
    }
    \endverbatim


    \table
        Property     | Description             | Required    | Default value
        velocityDirection | direction of flow  | no   | (0 0 0) / patchNormal
        phi            | flux field name     | no          | phi
    \endtable

See also
Foam::directionMixedFvPatchFields
Foam::basePressureFvPatchScalarField

SourceFiles
    pressureVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureVelocityFvPatchVectorFields_H
#define pressureVelocityFvPatchVectorFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/directionMixed/directionMixedFvPatchFields.H"
#include "fvMatrices/fvBlockMatrix/fvBlockMatrix.H"
#include "derivedFvPatchFields/basePressure/basePressureFvPatchScalarField.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class pressureVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class pressureVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
public:

    // Public data

        //- porous jump mode

        enum directionSpecification
        {
            userSpecified,
            zeroGradient,
            patchNormal,
            tangentialVelocity
        };

        static const NamedEnum<directionSpecification, 4>
            directionSpecificationNames_;

private:


    // Private data

        //- Name of the flux field
        word phiName_;

        //- Name of the p field - Needed only if p is coupled with U and
        //  name of p != p
        word pName_;

        //- directionMethod
        directionSpecification directionType_;

        //- velocity direction
        vector velDir_;

        //- velocity direction field
        vectorField velDirField_;

        //- Reference Frame Ptr
        coordinateFrame* coorFramePtr_;

        //- Print backflow information
        //  Requires processor sync - false by default
        const Switch backFlowReport_;

        //- Total area for backflow calculation
        scalar totalArea_;


    // Private Functions

        // is coupled
        bool isCoupledSolver() const;

        //- transform and set direction
        void setDirection();

        //- Set fraction based on normal
        void setTangentialFraction(const scalarField& phip);

        //- transform and set tangential velocity
        void transformAndSetTangentialVelocity();

        //- Implicit coefficients that are needed from grad(p) and laplacian(p)
        tmp<vectorField> implicitCoeff
        (
            const basePressureFvPatchScalarField&
        ) const;

        //- Calc and report backflow
        void reportBackFlow(const fvsPatchField<scalar>&) const;

public:

    //- Runtime type information
    TypeName("pressureVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureVelocityFvPatchVectorField onto a new patch
        pressureVelocityFvPatchVectorField
        (
            const pressureVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureVelocityFvPatchVectorField
        (
            const pressureVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureVelocityFvPatchVectorField
        (
            const pressureVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

            //- Return true: this patch field is altered by assignment
            //  Note: Assigment operator changes slightly results
            virtual bool assignable() const
            {
                return true;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);

        // Access

            //- Return the name of the velocity field
            virtual const word& phiName() const
            {
                return phiName_;
            }

        //- Operators specific to NS-eq.

            //- Momentum - gradp BC couplings (like total-pressure BC)
            virtual void addMomentumGradPCoupledBC
            (
                fvBlockMatrix<vector>&,
                const volScalarField&
            );

            //- Momentum - laplacian BC couplings (like total-pressure BC)
            virtual void addContinuityCoupledBC
            (
                BlockLduSystem<vector, scalar>&,
                const volScalarField&,
                const surfaceTensorField&
            ) const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Use weak boundary type for tangential velocity
        //  set value coefficient to zero for outflow
        // virtual tmp<Field<vector>> valueBoundaryCoeffs
        // (
        //     const tmp<scalarField>&
        // ) const;

        //- Manipulate matrix with given weights
        virtual void boundaryRelaxMatrix(fvBlockMatrix<vector>&) const;

        //- Write
        virtual void write(Ostream&) const;

    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
