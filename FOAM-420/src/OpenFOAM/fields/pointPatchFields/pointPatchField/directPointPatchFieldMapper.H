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
    (c) 2013 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::directPointPatchFieldMapper

Description
    direct pointPatchFieldMapper

\*---------------------------------------------------------------------------*/

#ifndef directPointPatchFieldMapper_H
#define directPointPatchFieldMapper_H

#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapper.H"
#include "fields/Fields/fieldMappers/directFieldMapper/directFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class directPointPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class directPointPatchFieldMapper
:
    public pointPatchFieldMapper,
    public directFieldMapper
{
    // Private Data

        //- Addressing from new back to old
        const labelUList& addressing_;

        //- Does map contain any unmapped values
        bool hasUnmapped_;


public:

    // Constructors

        //- Construct given addressing
        directPointPatchFieldMapper(const labelUList& addressing)
        :
            addressing_(addressing),
            hasUnmapped_(false)
        {
            if (addressing_.size() && min(addressing_) < 0)
            {
                hasUnmapped_ = true;
            }
        }


    //- Destructor
    virtual ~directPointPatchFieldMapper()
    {}


    // Member Functions

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override
        {
            return hasUnmapped_;
        }

        //- Return the direct map addressing
        virtual const labelUList& addressing() const override
        {
            return addressing_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
