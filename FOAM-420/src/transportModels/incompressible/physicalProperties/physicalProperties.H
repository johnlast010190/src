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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::physicalProperties

Description
    A container class for invariant incompressible physicalProperties
    Should be subsumed into universal thermo model

SourceFiles
    physicalProperties.C
    newphysicalProperties.C

\*---------------------------------------------------------------------------*/

#ifndef physicalProperties_H
#define physicalProperties_H

#include "db/dictionary/dictionary.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "memory/autoPtr/autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class physicalProperties Declaration
\*---------------------------------------------------------------------------*/

class physicalProperties
{

protected:

    // Protected data

        autoPtr<dimensionedScalar> rho_;
        autoPtr<dimensionedScalar> Cp_;
        autoPtr<dimensionedScalar> lambda_;
        autoPtr<dimensionedScalar> Prt_;
        autoPtr<dimensionedScalar> Pr_;

        dictionary physicalProperties_;

    // Private Member Functions

        //- Disallow copy construct
        physicalProperties(const physicalProperties&);

        //- Disallow default bitwise assignment
        void operator=(const physicalProperties&);


public:

        //- Construct from components
        physicalProperties
        (
            const dictionary& physicalProperties
        );


    // Destructor

        virtual ~physicalProperties()
        {}


    // Member Functions

        //- Return the physical properties dictionary
        const dictionary& physicalPropertiesDict() const
        {
            return physicalProperties_;
        }

        //- Return the density
        virtual const dimensionedScalar& rho() const;

        //- Return the constant pressure heat capacity
        virtual const dimensionedScalar& Cp() const;

        //- Return the thermal conductivity
        virtual const dimensionedScalar& lambda() const;

        //- Return the turbulent Prandtl number
        virtual const dimensionedScalar& Prt() const;

        //- Return the laminar Prandtl number
        virtual const dimensionedScalar& Pr() const;

        //- return autoPtr status
        inline virtual bool rhoValid() const
        {
            return rho_.valid();
        };

        //- return autoPtr status
        inline virtual bool CpValid() const
        {
            return Cp_.valid();
        }

        //- return autoPtr status
        inline virtual bool lambdaValid() const
        {
            return lambda_.valid();
        }

        //- return autoPtr status
        inline virtual bool PrtValid() const
        {
            return Prt_.valid();
        }

        //- return autoPtr status
        inline virtual bool PrValid() const
        {
            return Pr_.valid();
        }

        //- Correct the physical properties
        virtual void correct()
        {}

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& physicalProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
