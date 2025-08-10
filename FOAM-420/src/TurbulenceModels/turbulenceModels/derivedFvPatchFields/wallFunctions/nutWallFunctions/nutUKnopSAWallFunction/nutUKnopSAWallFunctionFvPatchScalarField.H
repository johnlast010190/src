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
    Foam::nutUKnopSAWallFunctionFvPatchScalarField

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
    nutUKnopSAWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUKnopSAWallFunctionFvPatchScalarField_H
#define nutUKnopSAWallFunctionFvPatchScalarField_H

#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUKnopSAWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUKnopSAWallFunctionFvPatchScalarField
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
        scalar CySA_;



protected:

    // Protected Member Functions

        //- simple y+ calc
        scalar calcYPlus
        (
            const scalar y,
            const scalar utau,
            const scalar nu
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
            const label N,
            const scalarField& utauInit
        ) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutUKnopSAWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUKnopSAWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUKnopSAWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUKnopSAWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUKnopSAWallFunctionFvPatchScalarField
        (
            const nutUKnopSAWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUKnopSAWallFunctionFvPatchScalarField
        (
            const nutUKnopSAWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUKnopSAWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUKnopSAWallFunctionFvPatchScalarField
        (
            const nutUKnopSAWallFunctionFvPatchScalarField&,
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
                new nutUKnopSAWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
