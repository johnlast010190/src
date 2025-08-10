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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "alembic/AlembicReader.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/gzstream/gzstream.h"
#include "include/OSspecific.H"


/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileFormats
{
    defineTypeNameAndDebug(AlembicReader, 0);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileFormats::AlembicReader::readFile
(
    const fileName& filename
)
{
// open the iarchive
    AbcF::IFactory factory;
    factory.setPolicy(Abc::ErrorHandler::kQuietNoopPolicy);
    AbcF::IFactory::CoreType coreType;
    Abc::IArchive archive = factory.getArchive(filename, coreType);
    //IArchive archive(Alembic::AbcCoreOgawa::ReadArchive(),argv[1]);

// check core type (must be Ogawa since HDF5 not compiled)
    checkCoreType(coreType);

// get parent object in hierarchy
    Abc::IObject obj=archive.getTop();

    if (debug)
    {
        // print file metadata
        getMetaData(archive);

        // print time series information
        printTimeSeriesData(archive);

        //traverse data hierarchy and print info
        traverseDataTree(obj, true);
    }

    // helper variables for assigning patch group names
    label groupID = 0;
    label maxGroupID = 0;

// recursively get all the mesh data from the archive
    for (uint32_t i=0; i<obj.getNumChildren(); i++)
    {
        Abc::IObject cobj = obj.getChild(i);

        const AbcA::MetaData& md = cobj.getMetaData();
        std::string schema = md.get("schema");

        if (AbcG::IXformSchema::matches(md))
        {
            AbcG::IXform xform(cobj, AbcG::kWrapExisting);
            DynamicList<AbcG::IPolyMesh> meshes;

            for (uint32_t j=0; j<cobj.getNumChildren(); j++)
            {
                Abc::IObject ccobj = cobj.getChild(j);
                const AbcA::MetaData& md2 = ccobj.getMetaData();
                std::string schema2 = md2.get("schema");

                if (AbcG::IPolyMeshSchema::matches(md2))
                {
                    meshes.append(AbcG::IPolyMesh(ccobj));

                    // add object name to group
                    word group = ccobj.getName();
                    HashTable<label>::const_iterator findGroup =
                        groupToPatch_.find(group);
                    if (findGroup != groupToPatch_.end())
                    {
                        groupID = findGroup();
                    }
                    else
                    {
                        groupID = maxGroupID;
                        groupToPatch_.insert(group, groupID);
                        maxGroupID++;
                    }
                }
            }
            meshes.shrink();

            if (meshes.size())
            {
                meshdata_.append(Tuple2(xform, meshes));
            }
        }
    }
    meshdata_.shrink();

    return true;
}


void Foam::fileFormats::AlembicReader::checkCoreType
(
    AbcF::IFactory::CoreType& coreType
)
{
    word coreName;
    if (coreType == AbcF::IFactory::kOgawa)
    {
        coreName = "Ogawa";
    }
    else if (coreType == AbcF::IFactory::kHDF5)
    {
        coreName = "HDF5";
        FatalErrorInFunction
            << "HDF5 format not supported, please use Ogawa"
            << exit(FatalError);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown .abc file format"
            << exit(FatalError);
    };
    if (debug)
    {
        Info<<"  core type : " << coreName <<"\n";
    }
}


void Foam::fileFormats::AlembicReader::getMetaData
(
    Abc::IArchive& archive
)
{
    std::string appName;
    std::string libraryVersionString;
    Alembic::Util::uint32_t libraryVersion;
    std::string whenWritten;
    std::string userDescription;
    std::string coreName;
    double dccFps;

    GetArchiveInfo
    (
        archive,
        appName,
        libraryVersionString,
        libraryVersion,
        whenWritten,
        userDescription,
        dccFps
    );

    if (appName != "")
    {
        Info<< "  file written by: " << appName << endl;
        Info<< "  using Alembic : " << libraryVersionString << endl;
        Info<< "  written on : " << whenWritten << endl;
        Info<< "  user description : " << userDescription << endl;
        if (dccFps > 0)
        {
            Info<< "  DCC FPS at write: " << dccFps << endl;
        }
    }
    else
    {
        Info<< "  (file doesn't have any ArchiveInfo)" << endl;
    }
}


void Foam::fileFormats::AlembicReader::printTimeSeriesData
(
    Abc::IArchive& archive
)
{
    uint32_t numTimes = archive.getNumTimeSamplings();

    for (uint32_t t=0; t<numTimes; t++)
    {
        AbcA::TimeSamplingPtr tsptr = archive.getTimeSampling(t);
        index_t maxSample =
            archive.getMaxNumSamplesForTimeSamplingIndex(t);

        Info<< "Time series " << t << " max Num Samples: "
            << maxSample << endl;

        for (index_t i=0; i<maxSample; i++)
        {
            Info<< "Index " << i << ": time = "
                << tsptr->getSampleTime(i) << endl;
        }
    }
}


void Foam::fileFormats::AlembicReader::printProperties
(
    Abc::ICompoundProperty iCProp,
    Abc::PropertyHeader header
)
{
    word ptype;
    AbcA::MetaData md;

    if (header.isCompound())
    {
        ptype = "CompoundProperty";
    }
    else if (header.isScalar())
    {
        ptype = "ScalarProperty";
    }
    else if (header.isArray())
    {
        ptype = "ArrayProperty";
    }
    else
    {
        ptype = "Unknown";
    }

    Info<< "property " << header.getName() << " type: " << ptype;

    if (header.isScalar())
    {
        Abc::IScalarProperty iProp(iCProp, header.getName());
        Info<< "[" << iProp.getNumSamples() << "]";
        md = iProp.getMetaData();
    }
    else if (header.isArray())
    {
        Abc::IArrayProperty iProp(iCProp, header.getName());
        Info<< "[" << iProp.getNumSamples() << "]";
        md = iProp.getMetaData();
    }

    //print metadata
    Info<< " {"<< md.serialize()<< "} " << endl;
}


void Foam::fileFormats::AlembicReader::printObjectInfo
(
    AbcG::IObject iObj
)
{
    AbcA::MetaData md = iObj.getMetaData();
    std::string schema = md.get("schema");
    Info<< "object " << iObj.getFullName() << " schema: " << schema;

    //print metadata
    Info<< " {"<< md.serialize()<< "} " << endl;
}


void Foam::fileFormats::AlembicReader::traverseCompoundProps
(
    Abc::ICompoundProperty iProp,
    bool first
)
{
    // header
    if (iProp.getNumProperties() > 0)
    {
        Info<< iProp.getObject().getFullName() << "/"
            << iProp.getName() << ":" << endl;
    }

    // children
    for (size_t c = 0; c < iProp.getNumProperties(); ++c)
    {
        printProperties(iProp, iProp.getPropertyHeader(c));
    }

    // visit children
    if (iProp.getNumProperties() > 0)
    {
        for (size_t p = 0; p < iProp.getNumProperties(); ++p)
        {
            Abc::PropertyHeader header = iProp.getPropertyHeader(p);
            if (header.isCompound())
            {
                traverseCompoundProps
                (
                    Abc::ICompoundProperty(iProp, header.getName())
                );
            }
        }
    }
}


void Foam::fileFormats::AlembicReader::traverseDataTree
(
    AbcG::IObject iObj,
    bool first
)
{
    Abc::ICompoundProperty props = iObj.getProperties();

    // header
    if (iObj.getNumChildren() > 0 || props.getNumProperties() > 0)
    {
        Info<< iObj.getFullName() << ":" << endl;
    }

    // children
    for (size_t c = 0; c < iObj.getNumChildren(); ++c)
    {
        printObjectInfo(iObj.getChild(c));
    }

    // properties
    for (size_t p = 0; p < props.getNumProperties(); ++p)
    {
        printProperties(props, props.getPropertyHeader(p));
    }

    // recursively visit property children
    if (props.getNumProperties() > 0)
    {
        for (size_t p = 0; p < props.getNumProperties(); ++p)
        {
            Abc::PropertyHeader header = props.getPropertyHeader(p);
            if (header.isCompound())
            {
                traverseCompoundProps
                (
                    Abc::ICompoundProperty(props, header.getName())
                );
            }
        }
    }

    // recursively visit object children
    if (iObj.getNumChildren() > 0)
    {
        for (size_t c = 0; c < iObj.getNumChildren(); ++c)
        {
            traverseDataTree(iObj.getChild(c));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::AlembicReader::AlembicReader
(
    const fileName& filename,
    bool storeBase
)
:
    points_(),
    groupToPatch_(),
    faces_(),
    meshdata_(),
    storeBasePoints_(storeBase),
    basePoints_(),
    mats_()
{
    readFile(filename);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileFormats::AlembicReader::~AlembicReader()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fileFormats::AlembicReader::clear()
{
    points_.clear();
    faces_.clear();
    basePoints_.clear();
    mats_.clear();
}


void Foam::fileFormats::AlembicReader::setTimeData(index_t timeIdx)
{
    // create sample selector from time index
    Abc::ISampleSelector iss(timeIdx);
    setTimeData(iss);
}


void Foam::fileFormats::AlembicReader::setTimeData
(
    Abc::chrono_t time,
    Abc::ISampleSelector::TimeIndexType idxType
)
{
    // create sample selector from time and index type
    Abc::ISampleSelector iss(time, idxType);
    setTimeData(iss);
}


void Foam::fileFormats::AlembicReader::setTimeData(Abc::ISampleSelector& iss)
{
    // clear existing mesh data
    clear();

    // offset for face indices
    label offset = 0;

    forAll(meshdata_, i)
    {
        // get transform data from IXform object
        AbcG::IXformSchema xschema = meshdata_[i].first().getSchema();
        AbcG::XformSample xform_samp;
        xschema.get(xform_samp, iss);
        Abc::M44d mat = xform_samp.getMatrix();

        DynamicList<AbcG::IPolyMesh>& meshes = meshdata_[i].second();
        forAll(meshes, j)
        {
            // get mesh data from IPolyMesh object
            AbcG::IPolyMeshSchema schema = meshes[j].getSchema();
            AbcG::IPolyMeshSchema::Sample mesh_samp;
            schema.get(mesh_samp, iss);

            label fiStartIdx = 0; // face indices start index
            uint32_t psize=mesh_samp.getPositions()->size();
            uint32_t fsize=mesh_samp.getFaceCounts()->size();

            // get patch group index from name
            label groupID = groupToPatch_.lookup
            (
                word(meshes[j].getName()), label(-1)
            );

            for (uint32_t k=0; k<psize; k++)
            {
                // get mesh points and transform
                V3f p=mesh_samp.getPositions()->get()[k];
                if (storeBasePoints_)
                {
                    basePoints_.append(point(p.x, p.y, p.z));
                    mats_.append(mat);
                }
                p=p*mat;
                points_.append(point(p.x, p.y, p.z));
            }
            for (uint32_t k=0; k<fsize; k++)
            {
                // given face count, create tri faces
                uint32_t fc=mesh_samp.getFaceCounts()->get()[k];
                // copied from readOBJ: simple triangulation algorithm for poly meshes
                List<label> verts(fc, -1);
                forAll(verts, vi)
                {
                    verts[vi] = mesh_samp.getFaceIndices()->get()[fiStartIdx+vi] + offset;
                }
                for (label fp = 1; fp < verts.size() - 1; fp++)
                {
                    label fp1 = verts.fcIndex(fp);
                    labelledTri tri(verts[0], verts[fp], verts[fp1], groupID);
                    faces_.append(tri);
                }
                fiStartIdx += fc;
            }
            offset += psize;
        }
    }
}



// ************************************************************************* //
