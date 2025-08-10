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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2009 OpenCFD Ltd.

Class
    Foam::turbulentMassFluxConcentrationFvPatchScalarField

Description
    Fixed mass flux boundary condition for concentration.

SourceFiles
    turbulentMassFluxConcentrationFvPatchScalarField.C

Example
    patchName
    {
       type turbulentMassFluxConcentration;
       massFluxType timemassflux; //or massflux for constant value
       timemassflux table ( ( 0.0 0.0 ) ( 4.999 0.0 ) ( 5.000 1.0 ) ( 10.0 1.0 ));
       value uniform 1;
       rho rhoEff;
    }


\*---------------------------------------------------------------------------*/

#ifndef turbulentMassFluxConcentrationFvPatchScalarFields_H
#define turbulentMassFluxConcentrationFvPatchScalarFields_H

#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fields/fvPatchFields/basic/fixedGradient/fixedGradientFvPatchFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class turbulentMassFluxConcentrationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentMassFluxConcentrationFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

    // Data types

        //- Enumeration listing the possible hest source input modes
        enum massFluxType
        {
            massFluxS,
            massFluxTV
        };

protected:

    //- mass flux [kg/m2/s]
    scalarField m_;

    //- TV mass flux [kg/m2/s]
    autoPtr<Function1<scalar>> mt_;

private:

    //- species name
    word speciesName_;

    //- density name
    word rhoName_;

    //- laminar mass diffusivity
    dimensionedScalar Dlam_;

    //- Heat source type names
    static const NamedEnum<massFluxType, 2> massFluxTypeNames_;

    //- mass flux type
    massFluxType massflux_;

public:

    //- Runtime type information
    TypeName("turbulentMassFluxConcentration");


    // Constructors

        //- Construct from patch and internal field
        turbulentMassFluxConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentMassFluxConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentMassFluxConcentrationFvPatchScalarField onto
        //  a new patch
        turbulentMassFluxConcentrationFvPatchScalarField
        (
            const turbulentMassFluxConcentrationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentMassFluxConcentrationFvPatchScalarField
        (
            const turbulentMassFluxConcentrationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentMassFluxConcentrationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentMassFluxConcentrationFvPatchScalarField
        (
            const turbulentMassFluxConcentrationFvPatchScalarField&,
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
                new turbulentMassFluxConcentrationFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // Access
            //non-const access to specie mass flux
            scalarField& massFlux()
            {
                return m_;
            }

        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
