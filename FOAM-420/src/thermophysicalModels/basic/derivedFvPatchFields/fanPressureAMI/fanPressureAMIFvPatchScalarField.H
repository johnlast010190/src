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
    (c) 2017 OpenCFD Ltd
    (c) 2022 Esi Ltd.

Class
    Foam::fanPressureAMIFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition can be applied on AMI boundaries to add pressure
    drop based on fan curve

    The switch nonDimensional can be used for a non-dimensional table. It needs
    inputs rpm and dm of the fan.

    The nonDimensional flux (phi*) for the table is calculate as :

    \verbatim
        phi* = 120*phi/(PI^3*dm*rpm)
        where:
            dm is the mean diameter.
            rpm is the RPM of the fan.
    \endverbatim

    The non-dimensional pressure (Psi*):

    \verbatim
        Psi* = 2 deltaP/(rho*(sqr(PI*omega*dm)))
        where:
            deltaP is the pressure drop
    \endverbatim

    The non-dimensional table should be given as Psi* = F(phi*).

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name | no | phi
        rho          | density name if required | no | rho
        fanCurve     | function1               | yes
        outOfBounds  | out of bounds handling  | yes         |
        directionMode    | direction definitions type | no | plane
        flowDirection    | flow direction if directionMode is plane | yes |
        origin    | centre of cylinder if directionMode is cylindrical | yes |
        radialFlowDirection   | in or out of cylinder if directionMode is cylindrical | yes |
        nonDimensional | uses non-dimensional table | no | false
        rpm          | fan rpm for non-dimensional table | no | 0.0
        dm           | mean diameter for non-dimensional table | no | 0.0
    \endtable


    -if plane: the user must provide flow direction
    -if cylindrical: the user must provide origin and radialFlowDirection. If
    origin is not provided, the origin will be estimated as the average of all
    the points of the patch, which for a cylinder will be located inside it.

    Example of the boundary condition specification:
    \verbatim
    AMIside1
    {
        type            fanPressureAMI;
        patchType       cyclicAMI;
        jump            uniform 0;
        fanCurve        table ( ( 0 0 ) ( 0.2 10 ) ( 0.5 5 ) );
        flowDirection   ( 1 0 0 );
        relax 0.5;
        value           uniform 0;
    }

    AMIside2
    {
        type            fanPressureAMI;
        patchType       cyclicAMI;
        jump            uniform 0;
        flowDirection   ( 1 0 0 );
        fanCurve        table ( ( 0 0 ) ( 0.2 10 ) ( 0.5 5 ) );
        relax 0.5;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::fanFvPatchField
    Foam::totalPressureFvPatchScalarField
    Foam::interpolationTable

SourceFiles
   fanPressureAMIFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fanPressureAMIFvPatchScalarField_H
#define fanPressureAMIFvPatchScalarField_H

#include "fields/fvPatchFields/derived/fixedJumpAMI/fixedJumpAMIFvPatchField.H"
#include "referenceFrames/referenceFrameFvPatch/referenceFrameFvPatch.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fanPressureAMIFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fanPressureAMIFvPatchScalarField
:
    public fixedJumpAMIFvPatchField<scalar>,
    public referenceFrameFvPatch<vector>
{
public:

    // Public data

        //- direction mode

        enum directionMode
        {
            plane,
            cylindrical
        };

        static const NamedEnum<directionMode, 2> directionModeNames_;

private:

    // Private data

        //- porous method
        directionMode directionMode_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Tabulated fan curve
        autoPtr<Function1<scalar>> fanCurve_;

        //- flowDirection
        vector flowDir_;

        //- origin -- needed only by the cylindrical-based flow direction
        vector origin_;

        //-radial flow direction. towards the centre or outside the cylinder
        // in or out are allowed
        word radialFlowDir_;

        //- helper int to reverse signs master or slave based on flow direction
        int reverseFlow_;

        //- Swtich for using non-dimensional curve
        Switch nonDimensional_;

        // Parameters for non-dimensional table

            //- Fan rpm
            scalar rpm_;

            //- Fan mean diameter
            scalar dm_;

        //- Static head calculator class
        staticHead staticHead_;

        //- Relaxation factor for updating the pressure
        scalar relax_;

        //- Min jump limiter
        scalar minJump_;

        //- Max jump limiter
        scalar maxJump_;

        //- Initialise reverseFlow
    void initializeFlowDirection();

public:

    //- Runtime type information
    TypeName("fanPressureAMI");


    // Constructors

        //- Construct from patch and internal field
        fanPressureAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanPressureAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fanPressureAMIFvPatchScalarField
        //  onto a new patch
        fanPressureAMIFvPatchScalarField
        (
            const fanPressureAMIFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fanPressureAMIFvPatchScalarField
        (
            const fanPressureAMIFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fanPressureAMIFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fanPressureAMIFvPatchScalarField
        (
            const fanPressureAMIFvPatchScalarField&,
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
                new fanPressureAMIFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        virtual const scalar& rpm() const
        {
            return rpm_;
        }

        virtual scalar& rpm()
        {
            return rpm_;
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
