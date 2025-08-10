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
    (c) 2019 Esi Ltd.

Class
    Foam::epsilonWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulence dissipation wall constraint
    for low- and high-Reynolds number turbulence models.

    The condition can be applied to wall boundaries for which it
    - calculates \c epsilon and \c G
    - specifies the near-wall epsilon value

    where

    \vartable
        epsilon | turblence dissipation field
        G       | turblence generation field
    \endvartable

    The low-Re correction is activated by setting the entry
    \c lowReCorrection to 'on'; in this mode the model switches between
    laminar and turbulent functions based on the laminar-to-turbulent y+ value
    derived from kappa and E.  When the \c lowReCorrection is inactive, the
    wall function operates in high-Re mode.

Usage
    \table
        Property        | Description              | Required  | Default value
        Cmu             | model coefficient        | no        | 0.09
        kappa           | Von Karman constant      | no        | 0.41
        E               | model coefficient        | no        | 9.8
        lowReCorrection | Low-Re correction active | no        | off
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            epsilonWallFunction;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField
    Foam::omegaWallFunctionFvPatchScalarField

SourceFiles
    epsilonWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonEBLagWallFunctionFvPatchScalarField_H
#define epsilonEBLagWallFunctionFvPatchScalarField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;
class ellipticBlendingKE;

/*---------------------------------------------------------------------------*\
           Class epsilonWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonEBLagWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

        //- Local copy of turbulence epsilon field
        scalarField epsilon_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- List of averaging corner weights
        List<List<scalar>> cornerWeights_;

        //- Lag version of Elliptic Blending
        Switch lag_;

    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;

        //- Set the master patch - master is responsible for updating all
        //  wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //  multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return non-const access to an epsilon patch
        virtual epsilonEBLagWallFunctionFvPatchScalarField& epsilonPatch
        (
            const label patchi
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const turbulenceModel& turbulence,
            scalarField& epsilon0
        );

        //- Calculate the epsilon
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& epsilon
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }

public:

    //- Runtime type information
    TypeName("epsilonEBLagWallFunction");


    // Constructors

        //- Construct from patch and internal field
        epsilonEBLagWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        epsilonEBLagWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  epsilonWallFunctionFvPatchScalarField
        //  onto a new patch
        epsilonEBLagWallFunctionFvPatchScalarField
        (
            const epsilonEBLagWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        epsilonEBLagWallFunctionFvPatchScalarField
        (
            const epsilonEBLagWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonEBLagWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        epsilonEBLagWallFunctionFvPatchScalarField
        (
            const epsilonEBLagWallFunctionFvPatchScalarField&,
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
                new epsilonEBLagWallFunctionFvPatchScalarField(*this, iF)
            );
        }

    //- Destructor
    virtual ~epsilonEBLagWallFunctionFvPatchScalarField()
    {}


    // Member functions

        // Access

            //- Return non-const access to the master's epsilon field
            scalarField& epsilon(bool init = false);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            virtual void updateWeightedCoeffs(const scalarField& weights);

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<scalar>& matrix,
                const scalarField& weights
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
