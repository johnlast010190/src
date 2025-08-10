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
    (c) 2010-2016 Esi Ltd.

Class

Description

SourceFiles
    filmTemperatureCoupledBaffleMixedFvPathceScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef filmTemperatureCoupledBaffleMixedFvPathceScalarField_H
#define filmTemperatureCoupledBaffleMixedFvPathceScalarField_H

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchFields.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/boundaryKappa/boundaryKappa.H"
#include "fields/fvPatchFields/derived/mappedField/mappedFieldFvPatchFields.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "primitives/bools/Switch/Switch.H"
#include "thermoSingleLayer/thermoSingleLayer.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
   Class filmTemperatureCoupledBaffleMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class filmTemperatureCoupledBaffleMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public mappedPatchBase,
    public boundaryKappa
{
public:

    typedef Foam::regionModels::surfaceFilmModels::thermoSingleLayer
        filmModelType;

private:

    // Private Data

        //- Dictionary Patch Refference
//        const dictionary& dict_;

        //- Name of film region
        const word filmRegionName_;

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Retrieve film model from the database
        const filmModelType& filmModel() const;

public:

    //- Runtime type information
    TypeName("filmTemperatureCoupledBaffleMixed");

    // Constructors

        //- Construct from patch and internal field
        filmTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  filmTemperatureCoupledBaffleMixedFvPatchScalarField onto a new patch
        filmTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const filmTemperatureCoupledBaffleMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new filmTemperatureCoupledBaffleMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        filmTemperatureCoupledBaffleMixedFvPatchScalarField
        (
            const filmTemperatureCoupledBaffleMixedFvPatchScalarField&,
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
                new filmTemperatureCoupledBaffleMixedFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

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
