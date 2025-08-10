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
    (c) 2019 Esi Ltd.

Class
    Foam::sampledTriSurfaceToMesh

Description
    A sampledSurface defined by a triSurface which 'cuts' the mesh using the
    cuttingTriSurface algorithm.

Note
    Does not actually cut until update() called.

SourceFiles
    sampledTriSurfaceToMesh.C

\*---------------------------------------------------------------------------*/

#ifndef sampledTriSurfaceToMesh_H
#define sampledTriSurfaceToMesh_H

#include "sampledSurface/sampledSurface/sampledSurface.H"
#include "surface/cuttingTriSurface/cuttingTriSurface.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sampledTriSurfaceToMesh Declaration
\*---------------------------------------------------------------------------*/

class sampledTriSurfaceToMesh
:
    public sampledSurface,
    public cuttingTriSurface
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& fvMesh_;

        //- Coordinate frame
        coordinateFrame* coorFramePtr_;

private:

    // Private data

        //- If restricted to zones, name of this zone or a regular expression
        keyType zoneKey_;

        //- Triangulated faces or keep faces as is
        const bool triangulate_;

        //- Whether to perform planar cut and filter based on distance
        const bool planarCut_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- Is the surface moving?
        bool moveSurface_;

        //- motionPoints
        mutable pointField motionPoints_;

    // Private Member Functions

        //- Sample field on faces
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


        template<class Type>
        tmp<Field<Type>>
        interpolateField(const interpolation<Type>&) const;


public:

    //- Runtime type information
    TypeName("sampledTriSurfaceToMesh");


    // Constructors

        //- Construct from components
        sampledTriSurfaceToMesh
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const keyType& zoneKey = word::null,
            const bool triangulate = true,
            const bool planarCut = false
        );

        //- Construct from dictionary
        sampledTriSurfaceToMesh
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledTriSurfaceToMesh();


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();

        //- Write down surface point to be able to restart
        //  if the motion is incremental
        virtual bool write();

        //- Points of surface
        virtual const pointField& points() const
        {
            return cuttingTriSurface::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return cuttingTriSurface::surfFaces();
        }

        //- Const access to per-face zone/region information
        //  Could instead return meshCells or cellZoneId of the meshCells.
        virtual const labelList& zoneIds() const
        {
            return Foam::emptyLabelList;
        }

        //- Face area magnitudes
        virtual const vectorField& Sf() const
        {
            return cuttingTriSurface::Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return cuttingTriSurface::magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return cuttingTriSurface::Cf();
        }


        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return cuttingTriSurface::meshCells();
        }

        //- Sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const;


        //- Sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const;

        //- Sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const;

        //- Sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const;

        //- Sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const;


        //- Interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const;


        //- Interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const;

        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurface/sampledTriSurfaceToMesh/sampledTriSurfaceToMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
