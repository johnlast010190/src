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

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/sensor/sensor/sensor.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::sensor<Type>::sensor
(
    const fvMesh& mesh,
    const dictionary& sensorDict
)
:
    regIOobject
    (
        IOobject
        (
            sensorDict.lookup("name"),
            mesh.thisDb().instance(),
            mesh
        )
    ),
    mesh_(mesh),
    dict_(sensorDict),
    writeSensorField_(dict_.lookupOrDefault<Switch>("writeSensorField", false))
{}


template<class Type>
Foam::sensor<Type>::sensor(const sensor<Type>& se)
:
    tmp<sensor<Type>>::refCount(),
    regIOobject
    (
        IOobject
        (
            se.dict_.lookup("name"),
            se.mesh_.thisDb().instance(),
            se.mesh_
        )
    ),
    mesh_(se.mesh_),
    dict_(se.dict_),
    writeSensorField_(dict_.lookupOrDefault<Switch>("writeSensorField", false))
{
    if (debug) {Info<<"call copy constructor for sensor " << dict_.lookup("name") << " and type " << dict_.lookup("type") << endl;}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::sensor<Type>::~sensor()
{
    if (debug)
    {
        Info<<"call destructor for sensor " << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::sensor<Type>::value
(
    const GeometricField<Type, fvPatchField, volMesh>& volField
) const
{
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::sensor<Type>::valueField() const
{
    NotImplemented;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
void Foam::sensor<Type>::writeField()
{
    if (writeSensorField_ == false)
    {
        return;
    }

    volScalarField sensorField
    (
        IOobject
        (
            sensorType() + "_" + fieldName(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("sensor", dimless, 0.0),
        calculatedFvPatchField<scalar>::typeName
    );

    sensorField.primitiveFieldRef() = this->valueField();

    sensorField.write();
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sensor<Type>& se
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const sensor<Type>&)"
    );

    os.writeKeyword("sensor");
    se.dict_.write(os);

    return os;
}


template<class Type>
bool Foam::sensor<Type>::writeData(Ostream& os) const
{
    os.writeKeyword("sensor");
    os << dict_;

    return true;
}


// ************************************************************************* //
