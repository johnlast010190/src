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
    (c) 2018-1018 Esi Ltd.

Class
    Foam::meshControl

Description
    Top level control of mesh generation process

SourceFiles
    meshControl.C

\*---------------------------------------------------------------------------*/

#ifndef meshControl_H
#define meshControl_H

#include "db/dictionary/dictionary.H"
#include "primitives/bools/Switch/Switch.H"
#include "primitives/enums/NamedEnum.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class meshControl Declaration
\*---------------------------------------------------------------------------*/

class meshControl
{
public:

        //- Which performance setting
        enum meshModeType
        {
            QUALITY,
            BALANCED,
            FAST,
            DRYRUN
        };

        static const NamedEnum<meshModeType, 4> meshModeTypeNames;

        //- Which algorithm setting
        enum meshAlgoType
        {
            STANDARD,
            DUAL,
            EXTRUDE,
            SHELL
        };

        static const NamedEnum<meshAlgoType, 4> meshAlgoTypeNames;

        //- Which algorithm setting
        enum blockMeshType
        {
            IMPORT,
            BLOCK,
            AUTO
        };

        static const NamedEnum<blockMeshType, 3> blockMeshTypeNames;

private:

    //- Whether to perform mesh refinement
    Switch refine_;

    //- Whether to perform mesh snap
    Switch snap_;

    //- Whether to add surface layers
    Switch layers_;

    //- Whether to perform topology changes
    Switch topoChanges_;

    //- Meshing performance setting
    meshModeType mode_;

    //- Meshing algorithm type
    meshAlgoType algo_;

    //- Starting block mesh type
    blockMeshType block_;

    //- Whether to use dual mesh algorithm
    Switch dual_;

    //- Whether to use OpenVDB for castellation
    Switch vdb_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshControl(const meshControl&);

        //- Disallow default bitwise assignment
        void operator=(const meshControl&);

public:

    //- Runtime type information
    TypeName("meshControl");

    // Constructors

        //- Construct from meshing dictionary
        meshControl(const dictionary& dict);

    // Destructor

    virtual ~meshControl(){};

    // Member Functions

        // Access

            //- Whether to perform mesh refinement
            bool refine() const
            {
                return refine_;
            }

            //- Whether to perform mesh refinement
            bool snap() const
            {
                return snap_;
            }

            //- Whether to add surface layers
            bool layers() const
            {
                return layers_;
            }

            //- Whether running with dual mesh algorithm
            bool dual() const
            {
                return dual_;
            }

            //- Whether to use OpenVDB for castellation
            bool vdb() const
            {
                return vdb_;
            }

            //- Whether to perform topo changes
            bool topoChanges() const
            {
                return topoChanges_;
            }

            //- Performance method
            meshModeType mode() const
            {
                return mode_;
            }

            //- Mesh algorithm type
            meshAlgoType algorithm() const
            {
                return algo_;
            }

            //- Starting blockMesh type
            blockMeshType block() const
            {
                return block_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
