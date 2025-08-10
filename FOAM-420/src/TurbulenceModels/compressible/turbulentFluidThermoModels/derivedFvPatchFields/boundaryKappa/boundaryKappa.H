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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::boundaryKappa

Description
    Common functions used in temperature coupled boundaries.

    The thermal conductivity \c kappa may be obtained by the following methods:
      - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
        defined by 'kappaName'
      - 'fluidThermo' : use fluidThermo and default
        compressible::turbulenceModel to calculate kappa
      - 'solidThermo' : use solidThermo kappa()

    \par Keywords provided by this class:
      \table
        Property     | Description                | Required    | Default value
        kappaMethod  | Thermal conductivity method        | yes |
        kappaName    | Name of thermal conductivity field | no  | none
        alphaAni     | Name of the non-isotropic alpha (thermal diffusivity: Kappa/Cp)    | no  | Anialpha
        phaseName    | Name of phase                      | no  | word::null
      \endtable

Usage
    \verbatim
    nonIsotropicWall
    {
        ...
        kappaMethod     solidThermo;
        kappaName       none;
        alphaAni        Anialpha;
        ...
    }
    \endverbatim

SourceFiles
    boundaryKappa.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryKappa_H
#define boundaryKappa_H

#include "fields/Fields/scalarField/scalarField.H"
#include "primitives/enums/NamedEnum.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class basicThermo;

/*---------------------------------------------------------------------------*\
                   Class boundaryKappa Declaration
\*---------------------------------------------------------------------------*/

class boundaryKappa
{
public:

    // Public enumerations

        //- Type of supplied Kappa
        enum KMethodType
        {
            mtFluidThermo,
            mtSolidThermo,
            mtLookup
        };


protected:

    // Protected data

        static const NamedEnum<KMethodType, 3> KMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;

        //- Object registry
        const objectRegistry& obr_;

        //- How to get K
        KMethodType method_;

        //- Name of thermal conductivity field (if looked up from database)
        word kappaName_;

        //- Name of the non-Isotropic alpha (default: Anialpha)
        word alphaAniName_;

        //- Name of group for lookup of thermo dictionaries (default: word::null)
        word phaseName_;

public:

    // Constructors

        //- Construct from patch and K name
        boundaryKappa
        (
            const objectRegistry& obr,
            const fvPatch& patch,
            const word& calculationMethod,
            const word& kappaName,
            const word& alphaAniName,
            const word& groupName=word::null
        );

        //- Construct from patch and dictionary
        boundaryKappa
        (
            const objectRegistry& obr,
            const fvPatch& patch,
            const dictionary& dict,
            const word& groupName=word::null
        );

         //- Construct from patch and boundaryKappa
        boundaryKappa
        (
            const objectRegistry& obr,
            const fvPatch& patch,
            const boundaryKappa& base
        );


    // Member functions

        //- Method to obtain K
        word KMethodName() const
        {
            return KMethodTypeNames_[method_];
        }

        //- Method to obtain K
        KMethodType KMethod() const
        {
            return method_;
        }

        //- Name of thermal conductivity field
        const word& kappaName() const
        {
            return kappaName_;
        }

        //- Return kappa or nf & kappa & nf in scalar / tensor case respectively
        tmp<scalarField> kappa() const;

        //- Return nf*alphaEff or nf&alphaAni*Cp/Cpv in scalar / tensor case
        //  respectively. (Note: alpha = kappa/Cp; alphaEff = kappa/Cpv)
        tmp<vectorField> nfAlphaEff() const;

        //- Return nf*kappa or nf&kappaAni in scalar / tensor case respectively
        tmp<vectorField> nfKappa() const;

        // Return reference to thermo object for the specified phase
        const basicThermo& thermo() const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
