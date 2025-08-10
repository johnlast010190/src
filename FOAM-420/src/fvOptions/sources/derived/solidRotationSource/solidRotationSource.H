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
    (c) 2018-2021 Esi Ltd.
    (c) 2012-2013 OpenFOAM Foundation

Class
    Foam::fv::solidRotationSource

Description

    Supplies rotational convective heat flux for solid equation, using GRF

    This source is only intended to be used in zones with rotational symmetry,
    allowing rotation of solids to be simulated without mesh movement.

    Minor deviations from rotational symmetry (i.e. boundary faces not being
    exactly orthogonal to the rotational flux) are taken care of at the expense
    of conservation, while larger deviations will result in unphysical results.

    \heading Source usage

    Example usage:
    \verbatim
    disc_MRF
    {
        type            solidRotationSource;
        active          true;
        referenceFrames
        {
            all
            {
                cellZone        disk;
                selectionMode   cellZone;
                referenceFrame  discFrame;
            }
        }
        fields          ( h );
    }
    \endverbatim

    Convection scheme should be specified as
    \verbatim
    divSchemes
    {
        div(phiOmega,h) Gauss upwind;
    }
    \endverbatim

    Do not use 'bounded' as this is accounted for explicitly.

    See also: GRFSource

SourceFiles
    solidRotationSource.C

\*---------------------------------------------------------------------------*/
#ifndef solidRotationSource_H
#define solidRotationSource_H

#include "sources/derived/GRFSource/GRFSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class solidRotationSource Declaration
\*---------------------------------------------------------------------------*/

class solidRotationSource
:
    public GRFSource
{
public:

    //- Runtime type information
    TypeName("solidRotationSource");


    // Constructors

        //- Construct from components
        solidRotationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construct
        solidRotationSource(const solidRotationSource&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const solidRotationSource&) = delete;


    //- Destructor
    //virtual ~solidRotationSource()
    //{}


    // Member Functions

        //- Return rotation induced flux field for the entire domain
        virtual tmp<surfaceScalarField> phiOmega() const override;

        // Add explicit and implicit contributions

            //- Scalar
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
