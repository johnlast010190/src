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
    (c) 2016 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "particle/particle.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::particle::propertyList_  = Foam::particle::propertyList();
Foam::string Foam::particle::propertyTypes_ = Foam::particle::propertyTypes();

std::size_t Foam::particle::sizeofPosition() const
{
    return reinterpret_cast<const char *>(&facei_) - reinterpret_cast<const char *>(&position_);
}

std::size_t Foam::particle::sizeofFields() const
{
    return reinterpret_cast<const char *>(this) + sizeof(particle) - reinterpret_cast<const char *>(&position_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle(const polyMesh& mesh, Istream& is, bool readFields)
:
    mesh_(mesh),
    position_(),
    celli_(-1),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(-1),
    tetPti_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(-1),
    noBasePt_(false)
{
    if (is.format() == IOstream::ASCII)
    {
        is  >> position_ >> celli_;

        if (readFields)
        {
            is  >> facei_
                >> stepFraction_
                >> tetFacei_
                >> tetPti_
                >> origProc_
                >> origId_;
        }
    }
    else
    {
        if (readFields)
        {
            is.read(reinterpret_cast<char*>(&position_), sizeofFields());
        }
        else
        {
            is.read(reinterpret_cast<char*>(&position_), sizeofPosition());
        }
    }

    // Check state of Istream
    is.check(FUNCTION_NAME);
}


void Foam::particle::writePosition(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position_ << token::SPACE << celli_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&position_), sizeofPosition());
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const particle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.position_
            << token::SPACE << p.celli_
            << token::SPACE << p.facei_
            << token::SPACE << p.stepFraction_
            << token::SPACE << p.tetFacei_
            << token::SPACE << p.tetPti_
            << token::SPACE << p.origProc_
            << token::SPACE << p.origId_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.position()),
            p.sizeofFields()
        );
    }

    return os;
}


// ************************************************************************* //
