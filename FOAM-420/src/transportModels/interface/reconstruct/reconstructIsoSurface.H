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
    Foam::reconstructIsoSurface

SourceFiles
    reconstruct.C

Author
    Rettenmaier adapting the code bases of Kunkelmann

Description
    Reconstructs the interface with an 0.5 alpha iso surface by interpolating
    cell centered alpha values to the cell points. At edges of the cell a
    0.5-alpha point is interpolated linearely if the alpha value from the two
    edge defining points is higher and lower then 0.5.
    Calculates the a distance field to the interface in a band of cells near
    the interface (if selected)

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Rettenmaier adapting the code bases of Kunkelmann (main developer).

    Method Development and Intellectual Property :
        Rettenmaier<rettenmaier@gsc.tu-darmstadt.de>
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


#ifndef reconstructIsoSurface_H
#define reconstructIsoSurface_H

#include "reconstruct/reconstruct.H"

#include "meshes/polyMesh/syncTools/syncTools.H"

#include "interpolation/volPointInterpolation/volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class reconstructIsoSurface Declaration
\*---------------------------------------------------------------------------*/

class reconstructIsoSurface
:
    public reconstruct
{

private:

    //- const reference to the mesh
    const fvMesh& mesh_;

    //- distance to interface
    //  distance > 0: denser fluid
    volScalarField distance_;

    //- points of the interface
    IOList<vector> interfacePoints_;        //TODO not unique with IsoSurface TODO

    //- centres of the polygons of the reconstructed interface
    volVectorField interfaceCenters_;

    //- maximum distance of the distance field to the interface
    scalar distanceThreshold_;

    //- switches the distance calculation on or off
    Switch calcSignedDistance_;

    //- interpolator for alpha from cells to points
    //autoPtr<alphaInterpolation> alphaInterpolator_;

    //- interpolated values of alpha at mesh points
    scalarField alphaP_;

    //- instance of distributeField class
    distributeField distributeField_;

    // Private Member Functions

        //- calculates nHatv interfaceDensity and interfacePoints
        void interfacePosition();

        //- transports nHatv away from the interface and calculates
        //  the distance of cellpoints to the interface
        void propagateNHatv();

public:

    // Runtime type name
        TypeName("reconstructIsoSurface");

    // Constructors

        //- Construct from components
        reconstructIsoSurface
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& transpProp,
            const List<bool>& isWallPatch,
            const volScalarField& isInterface
        );

    //- Destructor
    virtual ~reconstructIsoSurface(){}

    // Member Functions

        //- calculates interface normal nHatv and nHatfv
        void reconstructInterface();

    //- returns nHatf projected on the face normals
    //  OVERWRIDES nHatf() in reconstruct.H
    //  The via isoSurface reconstructed normal is not fit to be used
    //  for interfaceCompression, therefore the normal is obtained by
    //  gradAlpha
    const surfaceScalarField& nHatf()                                   //TODO calculate gradAlpha only once per timestep (see contactLine calc.)
    {
//        const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
//        surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
//        nHatf_ = gradAlphaf/(mag(gradAlphaf) + deltaN_) & alpha1_.mesh().Sf();

          nHatf_  = nHatfv_ & alpha1_.mesh().Sf();
        return nHatf_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
