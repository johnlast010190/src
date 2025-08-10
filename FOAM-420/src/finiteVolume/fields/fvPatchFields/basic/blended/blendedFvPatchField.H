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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

Class
    Foam::blendedFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a base class for 'blended' type boundary
    conditions, i.e. conditions that mix fixed value and patch-normal gradient
    conditions.

    The respective contributions from each is determined by a weight field:

        \f[
            x_p = w x_p + (1-w) \left(x_c + \frac{\nabla_\perp x}{\Delta}\right)
        \f]

    where
    \vartable
        x_p   | patch values
        x_c   | patch internal cell values
        w     | weight field
        \Delta| inverse distance from face centre to internal cell centre
        w     | weighting (0-1)
    \endvartable


Usage
    \table
        Property     | Description             | Required    | Default value
        valueFraction | weight field           | yes         |
        refValue     | fixed value             | yes         |
        refGrad      | patch normal gradient   | yes         |
    \endtable

Note
    This condition is not usually applied directly; instead, use a derived
    blended condition such as \c inletOutlet

See also
    Foam::inletOutletFvPatchField

SourceFiles
    blendedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef blendedFvPatchField_H
#define blendedFvPatchField_H

#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class blendedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blendedFvPatchField
:
    public fvPatchField<Type>
{
protected:

    // Protected data

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;

        //- Sensor used to determine valueFraction field
        //autoPtr<sensor<scalar>> bcSensor_;
        word sensorName_;

        //- functional relationship between input and output quantity
        autoPtr<Function1<scalar>> bcValueController_;

        //- switch already updated this time step!?
        bool updatedSwitch_;

        //- store time index at update switch time
        label timeIndex_=0;

        //- First boundary condition
        autoPtr<fvPatchField<Type>> boundaryOne_;

        //- Second boundary condition
        autoPtr<fvPatchField<Type>> boundaryTwo_;

        bool constructorCall_ = true;

    // Protected member functions

        //- Calculate valueFraction from sensor and Function1
        void calcSwitch();

public:

    //- Runtime type information
    TypeName("blended");


    // Constructors

        //- Construct from patch and internal field
        blendedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blendedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given blendedFvPatchField onto a new patch
        blendedFvPatchField
        (
            const blendedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        blendedFvPatchField
        (
            const blendedFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            if (debug)
            {
                Info<<"calling clone function for field : " << this->internalField().name() << endl;
            }
            return tmp<fvPatchField<Type>>
            (
                new blendedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blendedFvPatchField
        (
            const blendedFvPatchField<Type>&,
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
                new blendedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- for boundaries of blended type only
            virtual List<bool> fixesValues() const
            {
                if (this->size())
                {
                    List<bool> fixesValue(this->size(), false);
                    forAll(fixesValue, i)
                    {
                        if (valueFraction_[i] >= 0.5)
                        {
                            fixesValue[i] = boundaryOne_().fixesValue();
                        }
                        else
                        {
                            fixesValue[i] = boundaryTwo_().fixesValue();
                        }
                    }

                    return fixesValue;
                }

                return List<bool>(0, false);
            }

            //- Return true if this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }

            //- for boundaries of blended type only
            virtual List<bool> assignables() const
            {
                if (this->size())
                {
                    List<bool> assignable(this->size(), false);
                    forAll(assignable, i)
                    {
                        if (valueFraction_[i] >= 0.5)
                        {
                            assignable[i] = boundaryOne_().assignable();
                        }
                        else
                        {
                            assignable[i] = boundaryTwo_().assignable();
                        }
                    }

                    return assignable;
                }

                return List<bool>(0, false);
            }


        // Return defining fields

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }

            virtual const word& sensorName() const
            {
                return sensorName_;
            }

            virtual const autoPtr<Function1<scalar>> bcValueController() const
            {
                return bcValueController_;
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
                const fvPatchField<Type>&,
                const labelList&
            );

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMapGIB(const gibFvPatchFieldMapper&);


        // Evaluation functions

            //- Return the first BC
            const fvPatchField<Type>& boundaryOne() const
            {
                return boundaryOne_();
            }

            //- Return the second BC
            const fvPatchField<Type>& boundaryTwo() const
            {
                return boundaryTwo_();
            }

            //- Return the first BC ref
            fvPatchField<Type>& boundaryOne()
            {
                return boundaryOne_();
            }

            //- Return the second BC ref
            fvPatchField<Type>& boundaryTwo()
            {
                return boundaryTwo_();
            }

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
