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
    (c) 2017 Esi Ltd.
    (c) 2014 Hannes Kroeger

Class
    Foam::decayingTurbulenceFvPatchVectorField

Description
    direction: It is +1 or -1 depending on wheather the flow goes in positive or negative x-direction,respectively.
    LField: The integral length field.
    refField: The external velocity field.
    RField: The Reynolds stress tensor field.
    value: The value of the field; at the initial time it should be equal to the refField field.

    inlet
    {
    type decayingTurbulenceInflowGenerator;
     direction 1;
    LField nonuniform List<scalar> ...;
    refField nonuniform List<vector> ...;
    RField nonuniform List<symmTensor> ...;
    value nonuniform List<vector> ...;
    }

SourceFiles
    decayingTurbulenceFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef decayingTurbulenceFvPatchVectoR_H
#define decayingTurbulenceFvPatchVectoR_H

#if !defined( WIN32 ) && !defined( WIN64 )
#include "primitives/random/random11/random11.H"
#else
#include "primitives/random/Random/Random.H"
#endif

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "containers/LinkedLists/user/SLList.H"

#include "fields/fvPatchFields/derived/inflowGenerator/decayingVorton.H"

namespace Foam
{


class decayingTurbulenceFvPatchVectorField
: public fixedValueFvPatchField<vector>
{
protected:

    //- length scale field
    Field<scalar>          L_;

    //- mean velocity field
    Field<vector>          U_;

    //- mean stress field
    Field<symmTensor>      R_;

    //- time index to prevent multiple updates in single iteration
    label                  curTimeIndex_;

    //- discrete vortical structures
    SLList<decayingVorton> vortons_;

    //- Lund tensor based on mean stress
    Field<tensor>          Lund_;

    //- current stress tensor
    Field<symmTensor>      Rinst_;

    //- relaxation coefficient
    label                  ind_;

    //- inlet direction
    label                  direction_;

    //- vortex size cut-off
    scalar                 minVortonLength_;


    void initialiseFromDictionary(const dictionary&);

public:

    TypeName("decayingTurbulenceInflowGenerator");

    decayingTurbulenceFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const decayingTurbulenceFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const decayingTurbulenceFvPatchVectorField&
    );

    decayingTurbulenceFvPatchVectorField
    (
        const decayingTurbulenceFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchField<vector>> clone() const
    {
        return tmp<fvPatchField<vector>>(new decayingTurbulenceFvPatchVectorField(*this));
    }

    virtual tmp<fvPatchField<vector>> clone(const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchField<vector>>(new decayingTurbulenceFvPatchVectorField(*this, iF));
    }

    virtual void autoMap(const fvPatchFieldMapper&);

    virtual void rmap(const fvPatchField<vector>&, const labelList&);

    //- Map (and resize as needed) from self given a mapping object
    virtual void autoMapGIB(const gibFvPatchFieldMapper&);

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;

private:

    void doUpdate();

};
}

#endif
