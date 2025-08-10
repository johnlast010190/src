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
    (c) 2020-2024 Esi Ltd.

Class
    Foam::staticHead

InNamespace
    Foam

Description
    Add static head in total pressure based computations

Usage
    \table
        Property     | Description                | Required | Default value
        rho          | Density field name         | no       | rho
        hRefPoint    | Reference point where static head pressure = 0 | no | (see below)
        rhoRef       | Reference density          | (see below) |
        integrateStaticHead | Integrate rho*g*h to determine static pressure profile | no | off
        nNonOrthogonalCorrectors | Non-orthogonal correctors to use when integrating static head | no | 0
    \endtable

    The hydrostatic pressure head \int_hRef^h rho*g dy
    is computed and added to the computed pressure. Here hRef is supplied from
    the hRefPoint if specified. rho is the patch density field if
    integrateStaticHead is set, otherwise it is the constant value rhoRef if
    specified, and if not, the average density over the patch is
    used.

    If integrateStaticHead is set, hRefPoint must lie within the range of
    heights spanned by the patch (outside this range, there would be no
    density available to integrate). If unspecified, it defaults to the 'top'
    of the patch (minimum value of g & Cf over the patch).

    If integrateStaticHead is not set and hRefPoint is not specified, then it
    defaults to the hRef supplied by the solver if present, otherwise the
    origin.

    Finally, if the name of the presure field is p_rgh, rho*g*(h-hRef)
    is subtracted from the computed pressure to give the above static pressure
    distribution for physical pressure.


SourceFiles
    addStaticHead.C

\*---------------------------------------------------------------------------*/

#ifndef staticHead_H
#define staticHead_H

#include "fields/Fields/primitiveFields.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchField.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class staticHead Declaration
\*---------------------------------------------------------------------------*/

// staticHead class definition
class staticHead
{
    // Private data
    private:

        //- Integrated static pressure head - stored so that equation can
        //  be progressively converged
        scalarField pStatic_;

        //- Reference to the patch class that instantiated staticHead
        const fvPatchField<scalar>& fvPf_;

        //- Fixed density used to determine pressure head
        const scalar rhoRef_;

        //- Use average patch density for rhoRef
        const Switch useRhoAveAsRef_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        const word rhoName_;

        //- Whether the static pressure head is integrated from rho*g*h
        //  using the non-uniform patch rho
        const Switch integrateStaticHead_;

        //- Point used for pressure head reference height if specified
        autoPtr<point> hRefPoint_;

        //- Number of non-ortho corrrectors in static head calculation
        const label nNonOrthCorr_;

        //- Add static head pressure in total pressure based computations
        Switch addStaticHead_;

        //- Number of horizontal planes over which density is averaged,
        //  when integrating static head
        const label nAveragingPlanes_;

    protected:
        //- Return buoyant density on the patch from thermo library
        virtual tmp<scalarField> patchRho();


    public:
        // Constructor declarations

            //- staticHead constructor when constructed from patch
            staticHead
            (
                const fvPatchField<scalar>&
            );

            //- staticHead constructor when constructed from
            //  the internal field and a dictionary
            staticHead
            (
                const fvPatchField<scalar>&,
                const fvPatch&,
                const dictionary&
            );

            //- staticHead constructor when constructed by
            //  mapping onto a new patch
            staticHead
            (
                const fvPatchField<scalar>&,
                const staticHead&,
                const fvPatchFieldMapper&
            );

            //- staticHead constructor when the patchField is either
            //  constructed as a copy
            //  or
            //  contructed as a copy setting the internal field reference
            staticHead
            (
                const fvPatchField<scalar>&,
                const staticHead&
            );

        // Destructor
            virtual ~staticHead() {};

        // Member Functions

            //- Return static pressure head
            Foam::tmp<Foam::scalarField> computeStaticHead();

            //- Add static pressure head to field if active
            void addStaticHead(scalarField&);

            //- Is static head activated?
            const Switch active() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap
            (
                const staticHead& sh, const labelList& addr
            );

            //- Map (and resize as needed) from self given a mapping object
            void autoMapGIB
            (
                const gibFvPatchFieldMapper&
            );


        //- Write
        void write
        (
            Ostream&
        ) const;


        word getRhoName() const;

}; // end staticHead class declaration


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
