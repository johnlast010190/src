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
    (c) 2023 Esi Ltd.

Class
    Foam::volumeMassFractions

Description
    Base class to get volume/mass fractions.

SourceFiles
    volumeMassFractions.C

\*---------------------------------------------------------------------------*/

#ifndef volumeMassFractions_H
#define volumeMassFractions_H

#include "fields/volFields/volFields.H"
#include "containers/Lists/PtrList/PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class volumeMassFractions Declaration
\*---------------------------------------------------------------------------*/

class volumeMassFractions
:
    public regIOobject
{

public:

    //- Run time type information
    TypeName("volumeMassFractions");


    // Constructors

        //- Read fraction fields
        volumeMassFractions
        (
            const objectRegistry& obr,
            const word& phaseName
        );


    //- Destructor
    virtual ~volumeMassFractions()
    {}


    // Member functions

        //- Return true for active phases/species
        virtual bool active(label phaseiSpeciesi) const = 0;

        //- Return the bool list of active phases/species
        virtual const List<bool>& active() const = 0;

        //- Return passive phase index
        virtual label passiveIndex() const
        {
            return -1;
        }

        //- Return the mass/volume fraction fields
        virtual PtrList<volScalarField>& fractions() = 0;

        //- Return the const mass/volume fraction fields
        virtual const PtrList<volScalarField>& fractions() const = 0;

        //- Write for data
        virtual bool writeData(Ostream& os) const override
        {
            NotImplemented;
        }

        //- Write
        virtual bool write(const bool valid) const override
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
