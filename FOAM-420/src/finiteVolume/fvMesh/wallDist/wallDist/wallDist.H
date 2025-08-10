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
    (c) 2016 OpenCFD Ltd.
    (c) 2015-2016 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

Class
    Foam::wallDist

Description
    Interface to run-time selectable methods to calculate the distance-to-wall
    and normal-to-wall fields.

    Example of the wallDist specification in fvSchemes:
    \verbatim
        wallDist
        {
            method meshWave;

            // Optional entry enabling the calculation
            // of the normal-to-wall field
            nRequired false;

            // Optional entry delaying wall distance update to every n steps
            // Default is 1 (update every step)
            updateInterval 5;
        }
    \endverbatim

See also
    Foam::patchDistMethod::meshWave
    Foam::patchDistMethod::Poisson
    Foam::patchDistMethod::advectionDiffusion

SourceFiles
    wallDist.C

\*---------------------------------------------------------------------------*/

#ifndef wallDist_H
#define wallDist_H

#include "meshes/MeshObject/MeshObject.H"
#include "fvMesh/wallDist/patchDistMethods/patchDistMethod/patchDistMethod.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class wallDist Declaration
\*---------------------------------------------------------------------------*/

class wallDist
:
    public MeshObject<fvMesh, UpdateableMeshObject, wallDist>
{
    // Private data

        //- Set of patch IDs
        const labelHashSet patchIDs_;

        //- Name for the patch set, e.g. "wall"
        const word patchTypeName_;

        //- Construction dictionary
        const dictionary dict_;

        //- Run-time selected method to generate the distance-to-wall field
        mutable autoPtr<patchDistMethod> pdm_;

        //- Distance-to-wall field
        mutable volScalarField y_;

        //- Flag to indicate if the distance-to-wall field is required
        mutable bool nRequired_;

        //- Normal-to-wall field
        mutable tmp<volVectorField> n_;

        //- Update wall distance every updateInterval_ steps
        const label updateInterval_;

        //- Flag to indicate whether the wall distance requires updating
        bool requireUpdate_;


    // Private Member Functions

        //
        labelHashSet findPatchIDs(const fvMesh&) const;

        //- Read PDM method (set to meshWave if not present)
        const dictionary& setPDM
        (
            const fvMesh& mesh,
            const word& patchTypeName
        );

        //- Construct the normal-to-wall field as required
        void constructn() const;

        //- Disallow default bitwise copy construct
        wallDist(const wallDist&);

        //- Disallow default bitwise assignment
        void operator=(const wallDist&);


public:

    // Declare name of the class and its debug switch
    ClassName("wallDist");


    // Constructors

        //- Construct from mesh, patch IDs and optional patch type name
        wallDist
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs,
            const word& patchTypeName = "patch"
        );

        //- Construct from mesh, patch IDs, default patch distance method
        //  and optional patch type name
        wallDist
        (
            const fvMesh& mesh,
            const word& defaultPatchDistMethod,
            const labelHashSet& patchIDs,
            const word& patchTypeName = "patch"
        );

        //- Construct from mesh and optional patch type name
        wallDist
        (
            const fvMesh& mesh,
            const word& patchTypeName = "wall"
        );


    //- Destructor
    virtual ~wallDist();


    // Member Functions

        //- Return the patchIDs
        const labelHashSet& patchIDs() const
        {
            return patchIDs_;
        }

        //- Return reference to cached distance-to-wall field
        const volScalarField& y() const
        {
            return y_;
        }

        //- Return reference to cached normal-to-wall field
        const volVectorField& n() const;

        //- Update the y-field when the mesh moves
        virtual bool movePoints();

        //- Update the y-field when the mesh changes
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
