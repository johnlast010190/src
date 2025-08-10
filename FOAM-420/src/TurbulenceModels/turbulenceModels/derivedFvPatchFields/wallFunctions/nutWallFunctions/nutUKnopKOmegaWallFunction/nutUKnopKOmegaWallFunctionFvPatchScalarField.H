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
    (c) 2010-2016 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::nutUKnopKOmegaWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions for smooth walls, based on velocity,  using
    Knop's blended model consistent formulation for Spalart-Allmaras

    T. Knop; "Model consistent universal wall-functions for RANS turbulence
    modelling", BAIL International Conf. on Boundary and Interior Layers, 2006


    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            nutUKnopSAWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutUKnopKOmegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUKnopKOmegaWallFunctionFvPatchScalarField_H
#define nutUKnopKOmegaWallFunctionFvPatchScalarField_H

#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUKnopKOmegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUKnopKOmegaWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
private:
    // Private member data

        //- Reichardt's law coefficients
        scalar Cr1_;

        scalar Cr2_;

        scalar Cr3_;

        scalar Cr4_;

        scalar Br1_;

        scalar CyR_;


        //- Spalding law coefficients
        scalar Bsa1_;

        //- Knop S-A y+ Blending coefficient
        scalar CyKOmega_;

        //- Relaxation factor for Newton's method
        scalar newtonsRelax_;

protected:

    // Protected Member Functions

        //- simple y+ calc
        scalar calcYPlus
        (
            scalar y,
            scalar utau,
            scalar nu
        ) const;

        //- Reichardt Utau
        tmp<scalarField> uTauReichardt
        (
            const scalarField& magGradU,
            const scalarField& utauInit
        ) const;

        //- Spalding N, Utau
        tmp<scalarField> uTauSpaldingN
        (
            const scalarField& magGradU,
            label N,
            const scalarField& utauInit
        ) const;

        //- Calculate the turbulence viscosity
        tmp<scalarField> calcNut() const override;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const override;


public:

    //- Runtime type information
    TypeName("nutUKnopKOmegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUKnopKOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUKnopKOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUKnopKOmegaWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUKnopKOmegaWallFunctionFvPatchScalarField
        (
            const nutUKnopKOmegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUKnopKOmegaWallFunctionFvPatchScalarField
        (
            const nutUKnopKOmegaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const override
        {
            return tmp<fvPatchScalarField>
            (
                new nutUKnopKOmegaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUKnopKOmegaWallFunctionFvPatchScalarField
        (
            const nutUKnopKOmegaWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const override
        {
            return tmp<fvPatchScalarField>
            (
                new nutUKnopKOmegaWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            tmp<scalarField> yPlus() const override;


        // I-O

            //- Write
            void write(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
