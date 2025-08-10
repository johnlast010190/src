/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM® : Professional Open-source CFD
|   o   O   o    |
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------

License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM® <http://www.openfoam.org/>.

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
    © 2018 ESI Ltd.

Author
    2018. Daniel Deising (Esi Ltd.). All rights reserved.
    2018. Nikolaos Magoulas (Esi Ltd.). All rights reserved.

InClass
    Foam:sensor

Description
    Sensor class for sensor-Controlled boundary condition. Returns one value
    for a given field based on a specified condition (e.g. average, max,
    cellValue).

    \verbatim
        sensor
        {
            fieldName    <name>; // name of vol<Type>Field used as sensor
                                 // field must be registered!
            type         <typeName>; // available: cellValue/set
            <... additional entries for sensor types ...>
        }
    \endverbatim

SourceFiles
    sensor.C
    sensorNew.C

\*---------------------------------------------------------------------------*/

#ifndef sensor_H
#define sensor_H

#include "db/dictionary/dictionary.H"
#include "fields/Fields/Field/Field.H"
#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "meshes/MeshObject/MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Forward declarations
    class Time;

    // Forward declaration of friend functions and operators
    template<class Type> class sensor;
    template<class Type> Ostream& operator<<(Ostream&, const sensor<Type>&);

/*---------------------------------------------------------------------------*\
                         Class sensor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class sensor
:
    public tmp<sensor<Type>>::refCount,
    public regIOobject
{

protected:

    // Protected data

        //- Autopointer to mesh
        const fvMesh& mesh_;

        //- Dictionary
        const dictionary dict_;

        //- Switch for writing the sensor field
        Switch writeSensorField_;


private:

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const sensor<Type>&);


public:

    //- Runtime type information
    TypeName("sensor")

    //- Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        sensor,
        mesh,
        (
            const fvMesh& mesh,
            const dictionary& sensorDict
        ),
        (mesh, sensorDict)
    );


    // Constructor

        //- Construct from dictionary
        sensor
        (
            const fvMesh& mesh,
            const dictionary& sensorDict
        );

        //- Copy constructor
        sensor(const sensor<Type>& se);

        virtual autoPtr<sensor<Type>> clone() const
        {
            if (debug) {Info<<"call clone() for sensor " << dict_.lookup("name") << " and type " << dict_.lookup("type") << endl;}
            return autoPtr<sensor<Type>>(sensor<Type>::New(mesh_, dict_, dict_.lookup("name")));
        }

    //- Selector

        static autoPtr<sensor<Type>> New
        (
            const fvMesh& mesh,
            const dictionary& sensorDict,
            const word name = "none"
        );


    //- Destructor
    virtual ~sensor();


    // Member Functions

        // Access

            //- Return sensor field name
            word fieldName() const
            {
                return dict_.lookup("fieldName");
            }


            word sensorType() const
            {
                return dict_.lookup("type");
            }


        // Evaluation

            //- Return single value as a function of (scalar) independent
            //  variable
            virtual Type value
            (
                const GeometricField<Type,
                fvPatchField,
                volMesh>& volField
            ) const;

            //- Return a scalarField populated with the values of the sensor
            virtual tmp<Field<scalar>> valueField() const;


        // from meshObject

            //- Move points
            virtual bool movePoints(){return true;}

            //- Update mesh
            virtual void updateMesh(const mapPolyMesh& mpm){}


        // I/O

            //- Write the sensor field
            void writeField();

            //- Ostream Operator
            friend Ostream& operator<< <Type>
            (
                Ostream& os,
                const sensor<Type>& se
            );

            //- Write in dictionary format
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSensor(Type)                                                       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(sensor<Type>, 0);                      \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        sensor<Type>,                                                          \
        mesh                                                                   \
    );


#define makeSensorType(SS, Type)                                               \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(sensorTypes::SS<Type>, 0);             \
                                                                               \
    sensor<Type>::addmeshConstructorToTable<sensorTypes::SS<Type>>             \
        add##SS##Type##meshConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/sensor/sensor/sensor.C"
    #include "cfdTools/general/sensor/sensor/sensorNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
