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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::sensor<Type>> Foam::sensor<Type>::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word name
)
{
/*
    typedef sensor<Type> sensorType;

    if (mesh.thisDb().foundObject<sensorType>(name))
    {
        //const sensor<Type>& sensor(mesh.thisDb().objectRegistry::template lookupObject<sensor<Type>>(name));
        //mesh.thisDb().objectRegistry::template lookupObject<sensor<Type>>(name);
        //const sensor<Type> &sensor =
        //return autoPtr<sensor<Type>>(&sensor);
        //typedef GeometricField<Type, fvPatchField, volMesh> volField;
        //const volField& sensor(mesh.thisDb().lookupObject<volField>("test"));

        if (debug)
        {
            Info<< " return autoPtr to existing object " << name << " of type "
                << sensorType::typeName << endl;
        }

        const sensorType& sensorT(mesh.thisDb().lookupObject<sensorType>(name));

        return autoPtr<sensorType>(&const_cast<sensorType&>(sensorT));
    }
    else
    {*/
        const word type
        (
            dict.lookupOrDefault<word>("type", "none")
        );

        if (debug)
        {
            Info<< " construct new autoPtr " << name << " of type"
                << type << endl;
        }

        const auto ctor = ctorTableLookup("sensor type", meshConstructorTable_(), type);
        return ctor(mesh, dict);
    //}
}


// ************************************************************************* //
