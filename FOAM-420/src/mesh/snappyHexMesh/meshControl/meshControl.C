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
\*---------------------------------------------------------------------------*/

#include "meshControl/meshControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(meshControl, 0);

template<>
const char* Foam::NamedEnum
<
    Foam::meshControl::meshModeType,
    4
>::names[] =
{
    "quality",
    "balanced",
    "fast",
    "dryrun"
};


template<>
const char* Foam::NamedEnum
<
    Foam::meshControl::meshAlgoType,
    4
>::names[] =
{
    "standard",
    "dual",
    "extrude",
    "shell"
};

template<>
const char* Foam::NamedEnum
<
    Foam::meshControl::blockMeshType,
    3
>::names[] =
{
    "import",
    "block",
    "auto"
};

} // End namespace Foam
const Foam::NamedEnum<Foam::meshControl::meshModeType, 4>
    Foam::meshControl::meshModeTypeNames;

const Foam::NamedEnum<Foam::meshControl::meshAlgoType, 4>
    Foam::meshControl::meshAlgoTypeNames;

const Foam::NamedEnum<Foam::meshControl::blockMeshType, 3>
    Foam::meshControl::blockMeshTypeNames;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshControl::meshControl(const dictionary& dict)
:
    refine_(dict.lookupOrDefault<Switch>("castellatedMesh",true)),
    snap_(dict.lookupOrDefault<Switch>("snap",true)),
    layers_(dict.lookupOrDefault<Switch>("addLayers",true)),
    topoChanges_(dict.lookupOrDefault<Switch>("allowTopoChanges",true)),
    mode_(QUALITY),
    vdb_(false)
{
    if (dict.found("meshMode"))
    {
        mode_ = meshModeTypeNames.read(dict.lookup("meshMode"));
    }

    if (dict.found("meshAlgorithm"))
    {
        algo_ = meshAlgoTypeNames.read(dict.lookup("meshAlgorithm"));
        if (algo_ == meshControl::DUAL)
        {
            dual_ = true;
        }
    }
    else
    {
        if (dict.found("dualMesh"))
        {
            dual_  = readBool(dict.lookup("dualMesh"));
            if (dual_)
            {
                algo_ = meshControl::DUAL;
            }
            else
            {
                algo_ = meshControl::STANDARD;
            }
        }
        else
        {
            dual_ = false;
            algo_ = meshControl::STANDARD;
        }
    }

    Info<<"\nUsing meshing algorithm: "<< meshAlgoTypeNames[algo_] << nl <<endl;

    Info<<"Using meshing method: "<< meshModeTypeNames[mode_] << nl <<endl;

    if (dict.found("VDBrefinement"))
    {
        vdb_ = readBool(dict.lookup("VDBrefinement"));

        if (vdb_)
        {
            if (!dict.found("blockData"))
            {
                FatalErrorInFunction
                    << "You have selected VDBrefinement but "
                    << "blockData has not been provided."
                    << exit(FatalError);
            }

            if (!dict.found("VDBdomain"))
            {
                FatalErrorInFunction
                    << "You have selected VDBrefinement but "
                    << "VDBdomain has not been provided."
                    << exit(FatalError);
            }
            else if
            (
                !dict.subDict("VDBdomain").found("min")
             || !dict.subDict("VDBdomain").found("max")
            )
            {
                FatalErrorInFunction
                    << "You have selected VDBrefinement but "
                    << "VDBdomain.min or VDBdomain.max has not been specified."
                    << exit(FatalError);
            }

            const scalar featAngle
            (
                readScalar
                (
                   dict.subDict("castellatedMeshControls").lookup
                   (
                      "resolveFeatureAngle"
                   )
                )
            );

            if (featAngle < 0 || featAngle > 180)
            {
                FatalErrorInFunction
                    << "For VDBrefinement please specify"
                    << "\n      0 < resolveFeatureAngle < 180"
                    << exit(FatalError);
            }
        } // if vdb_
    } // if found VDBrefinement


    if (vdb_)
    {
        block_ = meshControl::IMPORT;
    }
    else if (dict.found("blockType"))
    {
        block_ = blockMeshTypeNames.read(dict.lookup("blockType"));
    }
    else
    {
        if (dict.found("autoBlockMesh"))
        {
            bool autoBlock  = readBool(dict.lookup("autoBlockMesh"));
            if (autoBlock)
            {
                block_ = meshControl::AUTO;
            }
            else
            {
                block_ = meshControl::IMPORT;
            }
        }
        else
        {
            block_ = meshControl::IMPORT;
        }
    }

    Info<<"Using block method: "<< blockMeshTypeNames[block_] << nl <<endl;

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
