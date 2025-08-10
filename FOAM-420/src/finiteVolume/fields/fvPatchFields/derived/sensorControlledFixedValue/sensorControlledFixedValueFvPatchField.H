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

Class
    Foam::sensorControlledFixedValueFvPatchField

Description
    This boundary condition sets a uniform fixed value based on
    a controller table (e.g. set flowRate based on Temperature)
    and a given sensor type. Currently available sensors are
    cellValue and set. cellValue returns value of vol<Type>Field
    for given cell index and set return mean/min/max of vol<Type>Field
    over a given cellSet or cellZone.

    Example: set velocity based on temperature sensor
    \verbatim
        inlet
        {
            type     sensorControlledFixedValue;
            sensor
            {
                fieldName    T;
                type         set;
                setType      cellZone;
                cellZone     tempSensor;
                sensorOp     mean;
            }
            Function  table
            (
                (275.5 (0.0 0.0 0.0));
                (295.5 (1.1 2.2 3.3));
                (315.5 (4.4 5.5 6.6));
            )
            isInlet     true;
            value       uniform (0 0 0);
        }
    \endverbatim

SourceFiles
    sensorControlledFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef sensorControlledFixedValueFvPatchField_H
#define sensorControlledFixedValueFvPatchField_H

#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchField.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class sensorControlledFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class sensorControlledFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
protected:

    // Private data

        //- sensor which triggers changes in bc
        // NOTE: currently hard-coded to be a scalar-valued quantity!
        word sensorName_;

        //- name of sensor field
        word sensorFieldName_; //const

        //- functional relationship between input and output quantity
        autoPtr<Function1<Type>> bcValueController_;

        //- is inlet or outlet bc
        bool isInlet_;

public:

   //- Runtime type information
   TypeName("sensorControlledFixedValue");


   // Constructors

        //- Construct from patch and internal field
        sensorControlledFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sensorControlledFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  sensorControlledFixedValueFvPatchField
        //  onto a new patch
        sensorControlledFixedValueFvPatchField
        (
            const sensorControlledFixedValueFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        sensorControlledFixedValueFvPatchField
        (
            const sensorControlledFixedValueFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new sensorControlledFixedValueFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sensorControlledFixedValueFvPatchField
        (
            const sensorControlledFixedValueFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new sensorControlledFixedValueFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- check the sign for inlet/outlet
        void checkSign(Type& controllerValue) const
        {}

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// template specialization for vectors: check sign based on inlet/outlet
template<>
void sensorControlledFixedValueFvPatchField<vector>::checkSign
(
    vector& controllerValue
) const
{
    // get average patch normal vector (magnitude unimportant)
    vector patchNormal( gSum(this->patch().Sf()) );

    // check sign of value for inlet/outlet
    if ((controllerValue & patchNormal) > 0 && isInlet_)
    {
        Info<< "inlet bc but vector points outwards \n"
            << "switching sign of bc value" << endl;
        controllerValue = -controllerValue;
    }
    if ((controllerValue & patchNormal) < 0 && !isInlet_)
    {
        Info<< "outlet bc but vector points inwards \n"
            << "switching sign of bc value" << endl;
        controllerValue = -controllerValue;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
