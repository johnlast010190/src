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

\*---------------------------------------------------------------------------*/

#include "parcels/Templates/CollidingParcel/CollidingParcel.H"
#include "db/IOstreams/IOstreams.H"
#include "db/IOobjects/IOField/IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::CollidingParcel<ParcelType>::propertyList_ =
    Foam::CollidingParcel<ParcelType>::propertyList();

template<class ParcelType>
std::size_t Foam::CollidingParcel<ParcelType>::sizeofFields() const
{
    return reinterpret_cast<const char *>(&collisionRecords_) - reinterpret_cast<const char *>(&f_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    f_(Zero),
    angularMomentum_(Zero),
    torque_(Zero),
    collisionRecords_()
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> f_;
            is >> angularMomentum_;
            is >> torque_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&f_), sizeofFields());
        }

        is >> collisionRecords_;
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::CollidingParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<vector> f(c.fieldIOobject("f", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, f);

    IOField<vector> angularMomentum
    (
        c.fieldIOobject("angularMomentum", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, angularMomentum);

    IOField<vector> torque
    (
        c.fieldIOobject("torque", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, torque);

    labelFieldCompactIOField collisionRecordsPairAccessed
    (
        c.fieldIOobject("collisionRecordsPairAccessed", IOobject::MUST_READ),
        valid
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairAccessed);

    labelFieldCompactIOField collisionRecordsPairOrigProcOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigProcOfOther",
            IOobject::MUST_READ
        ),
        valid
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairOrigProcOfOther);

    labelFieldCompactIOField collisionRecordsPairOrigIdOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigIdOfOther",
            IOobject::MUST_READ
        ),
        valid
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairOrigProcOfOther);

    pairDataFieldCompactIOField collisionRecordsPairData
    (
        c.fieldIOobject("collisionRecordsPairData", IOobject::MUST_READ),
        valid
    );
    c.checkFieldFieldIOobject(c, collisionRecordsPairData);

    labelFieldCompactIOField collisionRecordsWallAccessed
    (
        c.fieldIOobject("collisionRecordsWallAccessed", IOobject::MUST_READ),
        valid
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallAccessed);

    vectorFieldCompactIOField collisionRecordsWallPRel
    (
        c.fieldIOobject("collisionRecordsWallPRel", IOobject::MUST_READ),
        valid
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallPRel);

    wallDataFieldCompactIOField collisionRecordsWallData
    (
        c.fieldIOobject("collisionRecordsWallData", IOobject::MUST_READ),
        valid
    );
    c.checkFieldFieldIOobject(c, collisionRecordsWallData);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        CollidingParcel<ParcelType>& p = iter();

        p.f_ = f[i];
        p.angularMomentum_ = angularMomentum[i];
        p.torque_ = torque[i];

        p.collisionRecords_ = collisionRecordList
        (
            collisionRecordsPairAccessed[i],
            collisionRecordsPairOrigProcOfOther[i],
            collisionRecordsPairOrigIdOfOther[i],
            collisionRecordsPairData[i],
            collisionRecordsWallAccessed[i],
            collisionRecordsWallPRel[i],
            collisionRecordsWallData[i]
        );

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::CollidingParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<vector> f(c.fieldIOobject("f", IOobject::NO_READ), np);
    IOField<vector> angularMomentum
    (
        c.fieldIOobject("angularMomentum", IOobject::NO_READ),
        np
    );
    IOField<vector> torque(c.fieldIOobject("torque", IOobject::NO_READ), np);

    labelFieldCompactIOField collisionRecordsPairAccessed
    (
        c.fieldIOobject("collisionRecordsPairAccessed", IOobject::NO_READ),
        np
    );
    labelFieldCompactIOField collisionRecordsPairOrigProcOfOther
    (
        c.fieldIOobject
        (
            "collisionRecordsPairOrigProcOfOther",
            IOobject::NO_READ
        ),
        np
    );
    labelFieldCompactIOField collisionRecordsPairOrigIdOfOther
    (
        c.fieldIOobject("collisionRecordsPairOrigIdOfOther", IOobject::NO_READ),
        np
    );
    pairDataFieldCompactIOField collisionRecordsPairData
    (
        c.fieldIOobject("collisionRecordsPairData", IOobject::NO_READ),
        np
    );
    labelFieldCompactIOField collisionRecordsWallAccessed
    (
        c.fieldIOobject("collisionRecordsWallAccessed", IOobject::NO_READ),
        np
    );
    vectorFieldCompactIOField collisionRecordsWallPRel
    (
        c.fieldIOobject("collisionRecordsWallPRel", IOobject::NO_READ),
        np
    );
    wallDataFieldCompactIOField collisionRecordsWallData
    (
        c.fieldIOobject("collisionRecordsWallData", IOobject::NO_READ),
        np
    );

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const CollidingParcel<ParcelType>& p = iter();

        f[i] = p.f();
        angularMomentum[i] = p.angularMomentum();
        torque[i] = p.torque();

        collisionRecordsPairAccessed[i] = p.collisionRecords().pairAccessed();
        collisionRecordsPairOrigProcOfOther[i] =
            p.collisionRecords().pairOrigProcOfOther();
        collisionRecordsPairOrigIdOfOther[i] =
            p.collisionRecords().pairOrigIdOfOther();
        collisionRecordsPairData[i] = p.collisionRecords().pairData();
        collisionRecordsWallAccessed[i] = p.collisionRecords().wallAccessed();
        collisionRecordsWallPRel[i] = p.collisionRecords().wallPRel();
        collisionRecordsWallData[i] = p.collisionRecords().wallData();

        i++;
    }

    const bool valid = (np > 0);

    f.write(valid);
    angularMomentum.write(valid);
    torque.write(valid);

    collisionRecordsPairAccessed.write(valid);
    collisionRecordsPairOrigProcOfOther.write(valid);
    collisionRecordsPairOrigIdOfOther.write(valid);
    collisionRecordsPairData.write(valid);
    collisionRecordsWallAccessed.write(valid);
    collisionRecordsWallPRel.write(valid);
    collisionRecordsWallData.write(valid);
}


template<class ParcelType>
template<class CloudType>
void Foam::CollidingParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    label np = c.size();

    IOField<vector>& f(cloud::createIOField<vector>("f", np, obr));
    IOField<vector>& angularMomentum
    (
        cloud::createIOField<vector>("angularMomentum", np, obr)
    );
    IOField<vector>& torque(cloud::createIOField<vector>("torque", np, obr));

    label i = 0;
    forAllConstIter(typename CloudType, c, iter)
    {
        const CollidingParcel<ParcelType>& p = iter();

        f[i] = p.f();
        angularMomentum[i] = p.angularMomentum();
        torque[i] = p.torque();

        i++;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CollidingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.f_
            << token::SPACE << p.angularMomentum_
            << token::SPACE << p.torque_
            << token::SPACE << p.collisionRecords_;
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.f_),
            p.sizeofFields()
        );
        os  << p.collisionRecords();
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
