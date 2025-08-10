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
    (c) held by original author

Class
    Foam::reconstruct

SourceFiles
    reconstruct.C

Author
    Rettenmaier

Description
    Base clase for explicite and implicte reconstructions of the interface
    which holds the interface normal and the interface normal interpolated on
    the cell phases

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Rettenmaier (main developer).

    Method Development and Intellectual Property :
        Rettenmaier
      Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
      Daniel Deising <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
      Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef reconstruct_H
#define reconstruct_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/volFields/volFields.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

#include "finiteVolume/fvc/fvc.H"
#include "interface/distributeField/distributeField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class reconstruct Declaration
\*---------------------------------------------------------------------------*/

class reconstruct
{

protected:

    // Protected data

        word name_;

        const volScalarField& alpha1_;

        //- interface normal
        volVectorField nHatv_;

        //- interpolated normal vectors
        surfaceVectorField nHatfv_;

        //- nHatf_ = nHatfv_ & alpha1_.mesh().Sf()
        surfaceScalarField nHatf_;

        //- interfacial area density (m^2/m^3)
        volScalarField interfaceDensity_;

        //- indicates if cell is treaten as interface
        const volScalarField& isInterface_;

        //- indicates if a boundary patch is of type wall
        const List<bool>& isWallPatch_;

         //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

    //    //- points of the interface
    //    List<vector> interfacePoints_;        // TODO not unique with IsoSurface

    //    //- distance to interface
    //    //  distance > 0: denser fluid
    //    volScalarField distance_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        reconstruct(const reconstruct&);

        //- Disallow default bitwise assignment
        void operator=(const reconstruct&);

public:

    // Runtime type name
        TypeName("reconstruct");

    // Declare run-time selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            reconstruct,
            dictionary,
            (
                const word& name,
                const volScalarField& alpha,
                const dictionary& transpProp,
                const List<bool>& isWallPatch,
                const volScalarField& isInterface
            ),
            (name, alpha, transpProp, isWallPatch, isInterface)
        );

    // Selectors

        //- Select null constructed
        static autoPtr<reconstruct> New
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& transpProp,
            const List<bool>& isWallPatch,
            const volScalarField& isInterface
        );

    // Constructors

        //- Construct from components
        reconstruct
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& transpProp,
            const List<bool>& isWallPatch,
            const volScalarField& isInterface
        );

    //- Destructor
        virtual ~reconstruct(){}

    // Member Functions

        //- calculates interface normal nHatv and nHatfv
        virtual void reconstructInterface(){}

        // Access
        word name() const
        {
            return name_;
        }

        volVectorField& nHatv()
        {
            return nHatv_;
        }

        virtual const surfaceScalarField& nHatf()
        {
            nHatf_  = nHatfv_ & alpha1_.mesh().Sf(); //TODO best practice?
            return nHatf_;
        }

        surfaceVectorField& nHatfv()
        {
            return nHatfv_;
        }

        const volScalarField& interfaceDensity() const
        {
            return interfaceDensity_;
        }

        const List<bool>& isWallPatch() const
        {
            return isWallPatch_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
