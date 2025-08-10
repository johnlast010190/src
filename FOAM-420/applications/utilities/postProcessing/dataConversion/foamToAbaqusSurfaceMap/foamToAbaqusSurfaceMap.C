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
    (c) 1991-2007 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

Application
    foamToAbaqusSurfaceMap

Description
    Generate Abaqus pressure loads mapping Foam results (based on AMI).
    The target (Abaqus inp file) needs to be given as input.
    mpirun -np 4 foamToAbaqusSurfaceMap -targetSurface "PatchName.inp"
    -pMean -latestTime -parallel

    -flipLoad it can be used to flip the load sign in output
    -flipNorm Abaqus and Foam normals need to point in opposite directions
     to ensure that the mapping works. This flag can be used to revert
     the normals during the mapping process.

Example
foamToAbaqusSurfaceMap -parallel -latestTime -targetSurface abaqusSurfMesh.inp
-scale "(0.001 0.001 0.001)" -translate "(0 0 400)"

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "AMIInterpolation/AMIInterpolation/AMIPatchToPatchInterpolation.H"
#include "primitives/strings/stringOps/stringOps.H"
#include "fields/volFields/volFields.H"
#include "db/IOobjectList/IOobjectList.H"
#include "knn.h"
#include <map>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// extract data for each column of a line -
static std::string readINPToken
(
    const string& line,
    const size_t& width,
    size_t& index
)
{
    size_t indexStart, indexEnd;

    indexStart = index;

    indexEnd = line.find(',', indexStart);

    index = indexEnd + 1;

    if (indexEnd == std::string::npos)
    {
    indexEnd = line.size();
        index = indexEnd;
    }

    return line.substr(indexStart, indexEnd - indexStart);
}

static scalar parseINPCoord(const string& s)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign))());
        scalar exponent = readScalar(IStringStream(s.substr(expSign+1))());

        if (s[expSign] == '-')
        {
            exponent = -exponent;
        }
        return mantissa*pow(10, exponent);
    }
    else
    {
        return readScalar(IStringStream(s)());
    }
}

namespace Foam
{


void constructKNNdata
(
    const primitivePatch& sourcePatch,
    std::map<word,KNN>& bndMaps,
    std::map<word,std::vector<label>>& faceIdMaps
)
{
    std::map<word,DynamicList<point>> bsource_xyz_;
    word sourceDummyName = "source";
    if (Pstream::parRun())
    {
        List<vectorField> posLoc(Pstream::nProcs());
        label myid=Pstream::myProcNo();

        DynamicList<point> bxyz;
        const vectorField& pcf = sourcePatch.faceCentres();

        forAll(pcf, fI)
        {
            const vector& pt= pcf[fI];
            bxyz.append(pt);
        }

        posLoc[myid]=bxyz;
        Pstream::allGatherList(posLoc);
        vectorField posGlob=ListListOps::combine<vectorField>(posLoc,accessOp<vectorField>());
        bsource_xyz_[sourceDummyName]=posGlob;
    }

    if (Pstream::parRun())
    {
        std::map<word,DynamicList<point>>::iterator it;

        for (it=bsource_xyz_.begin();it!=bsource_xyz_.end();it++)
        {
            word bname=it->first;
            std::vector<std::vector<scalar>> tData;
            const DynamicList<point> &bflist=it->second;
            for (label j=0;j<bflist.size();j++)
            {
                const point &pt=bflist[j];
                std::vector<scalar> pts(3);
                pts[0]=pt.x();
                pts[1]=pt.y();
                pts[2]=pt.z();
                tData.push_back(pts);
                faceIdMaps[bname].push_back(j);
            }

            bndMaps[bname].build(tData);
        }

        Pstream::barrier();
    }
    else
    {
        const vectorField& pcf = sourcePatch.faceCentres();
        std::vector<std::vector<scalar>> tData;

        forAll(pcf, fI)
        {
            const vector& pt= pcf[fI];
            std::vector<scalar> pts(3);
            pts[0]=pt.x();
            pts[1]=pt.y();
            pts[2]=pt.z();
            tData.push_back(pts);
            faceIdMaps[sourceDummyName].push_back(fI);
        }
        bndMaps[sourceDummyName].build(tData);
    }
}


void readINP
(
    const fvMesh& mesh,
    word surfName,
    pointField& targetPoints,
    faceList& targetfaces,
    labelList& eid
)
{

    DynamicList<point> dPoints;
    DynamicList<label> dIndices;
    DynamicList<face> dFaces;
    DynamicList<label> dEid;
    DynamicList<label> dPid;

    string status ="NULL";
    string eltype ="NOEL";
    label groupId = -1;

    // read done from the master node
    if (Pstream::master())
    {
    word surPath;
    if (Pstream::parRun())
    {
        surPath = mesh.time().path()/".."/"constant"/"triSurface"/surfName;
    }
    else
    {
        surPath = mesh.time().path()/"constant"/"triSurface"/surfName;
    }


    IFstream is(surPath);

    if (!is.good())
    {
        FatalErrorInFunction
        << "Cannot read file " << surPath
        << exit(FatalError);
    }

    while (is.good())
    {
        size_t linei = 0;
        string line;
        is.getLine(line);

         // skip empty lines and comments **
        if (line.empty() || (line[0] == '*'&& line[1] == '*'))
        {
        continue;
        }

        // skip all the *.inp section not related
        // to the grid topology
        if (status == "ELEMENT")
        {
            //this breaks at the end of the ELEMENT section
            //if it finds an uncommented line with a
            //keyword different from *ELEMENT
            if ((line[0] == '*')&&
                !(stringOps::upper(line).find("*ELEMENT")!=std::string::npos))
            {
                Info<<"End of the mesh section of the *.inp file."<<endl;
                Info<<"Stop reading *.inp file at line number: "<<is.lineNumber()<<"."<<endl;
                Info<<"Line: "<<line<<endl;
                break;
            }
        }

        //status NODE or ELEMENT
        if (stringOps::upper(line).find("*NODE")!=std::string::npos)
        {
        status ="NODE";
        continue;
        }
        else if (stringOps::upper(line).find("*ELEMENT")!=std::string::npos)
        {
        status ="ELEMENT";

        readINPToken(line, 8, linei);
        readINPToken(line, 8, linei);
        string patchstring = readINPToken(line, 8, linei);

        if (stringOps::upper(patchstring).find("ELSET=")!=std::string::npos)
        {
            // new patch (ELSET) found
            groupId++;
            std::size_t elsetpos = stringOps::upper(patchstring).find("ELSET=");
            string patchname = patchstring.substr(elsetpos+6);
            Info<< " name of the patch "<<groupId<<" is "<<patchname<<endl;
        }

        // find element type now, S3 and S4 supported
        if (stringOps::upper(line).find("S3", 14)!=std::string::npos)
        {
            eltype ="S3";
        }
        else if (stringOps::upper(line).find("S4", 14)!=std::string::npos)
        {
            eltype ="S4";
        }
        else
        {
            Info<< "Error: unknown element type found "<<nl
            << "at line: "<<is.lineNumber()<<endl;
        }
        continue;
        }

        if (status == "NODE")
        {
        linei = 0;
        label index =
            readLabel(IStringStream(readINPToken(line, 8, linei))());
        scalar x = parseINPCoord(readINPToken(line, 8, linei));
        scalar y = parseINPCoord(readINPToken(line, 8, linei));
        scalar z = parseINPCoord(readINPToken(line, 8, linei));

        // Info<<" Node index "<<index<<" coord  "<<x<<", "<<y<<", "<<z<<endl;

        dIndices.append(index);
        dPoints.append(point(x, y, z));
        }
        else if (status == "ELEMENT") //S4 case
        {
        linei = 0;
        if (eltype=="S4")
        {
            face newFace;
            newFace.setSize(4);
            label id   = readLabel(IStringStream(readINPToken(line, 8, linei))());
            newFace[0] = readLabel(IStringStream(readINPToken(line, 8, linei))());
            newFace[1] = readLabel(IStringStream(readINPToken(line, 8, linei))());
            newFace[2] = readLabel(IStringStream(readINPToken(line, 8, linei))());
            newFace[3] = readLabel(IStringStream(readINPToken(line, 8, linei))());

            dFaces.append(newFace);
            dPid.append(groupId);
            dEid.append(id);

            // Info<<" Face Index "<<id<<" cells  "<<newFace[0]<<", "
            //     <<newFace[1]<<", "
            //     <<newFace[2]<<", "
            //     <<newFace[3]<<endl;

        }
        else if (eltype=="S3")
        {
            face newFace;
            newFace.setSize(3);
            label id   = readLabel(IStringStream(readINPToken(line, 8, linei))());
            newFace[0] = readLabel(IStringStream(readINPToken(line, 8, linei))());
            newFace[1] = readLabel(IStringStream(readINPToken(line, 8, linei))());
            newFace[2] = readLabel(IStringStream(readINPToken(line, 8, linei))());

            dFaces.append(newFace);
            dPid.append(groupId);
            dEid.append(id);

            // Info<<" Face Index "<<id<<" cells  "<<newFace[0]<<", "
            //     <<newFace[1]<<", "
            //     <<newFace[2]<<endl;

        }

        }

    }
    } //end read master

    dPoints.shrink();
    dIndices.shrink();
    dFaces.shrink();
    dEid.shrink();
    dPid.shrink();

    // Build inverse mapping (index to point)
    Map<label> indexToPoint(2*dIndices.size());
    forAll(dIndices, i)
    {
    indexToPoint.insert(dIndices[i], i);
    }

    forAll(dFaces, i)
    {
    face& f = dFaces[i];
    forAll(f,fp)
    {
        f[fp] = indexToPoint[f[fp]];
    }
    }

    // Transfer DynamicLists to straight ones.
    pointField allPoints(dPoints.size());
    allPoints.transfer(dPoints);
    dPoints.clear();
    targetfaces = dFaces;
    targetPoints = allPoints;
    eid = dEid;


}

HashTable<wordHashSet> candidateObjects
(
    const IOobjectList& objects,
    const wordHashSet& supportedTypes,
    const bool specifiedFields,
    const wordHashSet& selectedFields
)
{
    // Special case = no fields
    if (specifiedFields && selectedFields.empty())
    {
        return HashTable<wordHashSet>();
    }

    HashTable<wordHashSet> usable = objects.classes();

    // Limited to types that we explicitly handle
    usable.retain(supportedTypes);

    // If specified, further limit to selected fields
    if (specifiedFields)
    {
        forAllIters(usable, iter)
        {
            iter.object().retain(selectedFields);
        }

        // Prune entries without any fields
        usable.filterValues
        (
            [](const wordHashSet& vals){ return !vals.empty(); }
        );
    }

    return usable;
}


template<class GeoField>
label readFields
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    PtrList<GeoField>& fields
)
{
    label nFields = 0;

    // Available fields of type GeomField
    const wordList fieldNames = objects.sortedNames(GeoField::typeName);

    fields.setSize(fieldNames.size());

    // Construct the fields
    for (const word& fieldName : fieldNames)
    {
        if (selectedFields.empty() || selectedFields.found(fieldName))
        {
            fields.set
            (
                nFields++,
                new GeoField(*(objects[fieldName]), mesh)
            );
        }
    }

    fields.setSize(nFields);

    return nFields;
}


template<class GeoField>
void writeFieldsAMI
(
    const primitivePatch& sP,
    const primitivePatch& tP,
    const bool& isSourceSurfacesPresent,
    const word& sourcePatchName,
    const scalar& loadDir,
    autoPtr<AMIPatchToPatchInterpolation> amiInterpolatePtr,
    const labelList& eid,
    PtrList<GeoField>& fields,
    const fileName& abaqusDir,
    simpleVTKWriter& sourceVTK,
    simpleVTKWriter& targetVTK
)
{
    forAll(fields, i)
    {
        writeFieldAMI
        (
            sP,
            tP,
            isSourceSurfacesPresent,
            sourcePatchName,
            loadDir,
            amiInterpolatePtr(),
            eid,
            fields[i],
            abaqusDir,
            sourceVTK,
            targetVTK
        );
    }
}

template<class GeoField>
void writeFieldsKNN
(
    const primitivePatch& sP,
    const primitivePatch& tP,
    const bool& isSourceSurfacesPresent,
    const word& sourcePatchName,
    const scalar& loadDir,
    std::map<word,KNN>& bndMaps,
    std::map<word,std::vector<label>>& faceIdMaps,
    const labelList& eid,
    PtrList<GeoField>& fields,
    const fileName& abaqusDir,
    simpleVTKWriter& sourceVTK,
    simpleVTKWriter& targetVTK
)
{
    forAll(fields, i)
    {
        writeFieldKNN
        (
            sP,
            tP,
            isSourceSurfacesPresent,
            sourcePatchName,
            loadDir,
            bndMaps,
            faceIdMaps,
            eid,
            fields[i],
            abaqusDir,
            sourceVTK,
            targetVTK
        );
    }
}


template<class GeoField>
void writeFieldAMI
(
    const primitivePatch& sP,
    const primitivePatch& tP,
    const bool& isSourceSurfacesPresent,
    const word& sourcePatchName,
    const scalar& loadDir,
    const AMIPatchToPatchInterpolation& amiInterpolate,
    const labelList& eid,
    GeoField& field,
    const fileName& abaqusDir,
    simpleVTKWriter& sourceVTK,
    simpleVTKWriter& targetVTK
)
{
    typedef typename GeoField::value_type GeoType;

    convertPressureToSIUnits(field); //pascal

    scalar flipLoad = 1;

    if (field.dimensions() == dimPressure)
    {
        // revert load only at pressure field (if needed)
        flipLoad = loadDir;
    }


    Field<GeoType> spField(sP.size(), pTraits<GeoType>::zero);
    label ifaces = 0;
    if (isSourceSurfacesPresent)
    {
        label patchID(field.mesh().boundaryMesh().findPatchID(sourcePatchName));
        const Field<GeoType>& pf = field.boundaryField()[patchID];
        forAll(pf, facei)
        {
            spField[ifaces] = flipLoad*pf[facei];
            ifaces++;
        }
    }
    else
    {
        forAll(field.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(field.mesh().boundary()[patchi]))
            {
                const Field<GeoType>& pf = field.boundaryField()[patchi];
                forAll(pf, facei)
                {
                    spField[ifaces] = flipLoad*pf[facei];
                    ifaces++;
                }
            }
        }
    }

    Field<GeoType> tpField( amiInterpolate.interpolateToTarget(spField) );

    // write abaqus inp file
    if (Pstream::master())
    {
        OFstream abaqusos
        (
            abaqusDir/
            "abaqusSurfaceData_" + word(field.name()) + "_" +
            field.mesh().time().timeName() + ".inp"
        );

        abaqusos << "**FOAM to ABAQUS load File" << endl;
        abaqusos << "**" << field.name() << endl;
        //abaqusos << "**"<< field.name() << " Units: "
        //         << field.dimensions() <<endl;
        abaqusos << "**Pressure Units: " << field.dimensions() <<endl;

        abaqusos<<"*DLOAD"<<endl;

        forAll(tP,faceID)
        {
            abaqusos
                << eid[faceID] << "," << field.name() << ","
                << tpField[faceID] << endl;
        }
    }

    sourceVTK.addFaceData(field.name(), spField);
    targetVTK.addFaceData(field.name(), tpField);
}


template<class GeoField>
void writeFieldKNN
(
    const primitivePatch& sP,
    const primitivePatch& tP,
    const bool& isSourceSurfacesPresent,
    const word& sourcePatchName,
    const scalar& loadDir,
    std::map<word,KNN>& bndMaps,
    std::map<word,std::vector<label>>& faceIdMaps,
    const labelList& eid,
    GeoField& field,
    const fileName& abaqusDir,
    simpleVTKWriter& sourceVTK,
    simpleVTKWriter& targetVTK
)
{
    typedef typename GeoField::value_type GeoType;

    convertPressureToSIUnits(field); //pascal

    scalar flipLoad = 1;

    if (field.dimensions() == dimPressure)
    {
        // revert load only at pressure field (if needed)
        flipLoad = loadDir;
    }


    Field<GeoType> spField(sP.size(), pTraits<GeoType>::zero);
    label ifaces = 0;
    if (isSourceSurfacesPresent)
    {
        label patchID(field.mesh().boundaryMesh().findPatchID(sourcePatchName));
        const Field<GeoType>& pf = field.boundaryField()[patchID];
        forAll(pf, facei)
        {
            spField[ifaces] = flipLoad*pf[facei];
            ifaces++;
        }
    }
    else
    {
        forAll(field.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(field.mesh().boundary()[patchi]))
            {
                const Field<GeoType>& pf = field.boundaryField()[patchi];
                forAll(pf, facei)
                {
                    spField[ifaces] = flipLoad*pf[facei];
                    ifaces++;
                }
            }
        }
    }

    word patchName = "source";

    std::map<word,std::vector<label>> ublist;
    const vectorField& pcf = tP.faceCentres();
    forAll(pcf, facei)
    {
        const point& pt= pcf[facei];
        std::vector<scalar> pts(3);
        pts[0]=pt.x();
        pts[1]=pt.y();
        pts[2]=pt.z();
        std::vector<label> nb;
        std::vector<scalar> dist;
        bndMaps[patchName].search
        (
            pts,
            nb,
            dist
        );
        label ib=nb[0];
        label jc=faceIdMaps[patchName][ib];
        ublist[patchName].push_back(jc);
    }

    Pstream::barrier();

    std::map<word,DynamicList<GeoType>> GeoTypeField;
    if (Pstream::parRun())
    {
        List<Field<GeoType>> posLoc(Pstream::nProcs());
        label myid=Pstream::myProcNo();

        DynamicList<GeoType> bxyz;

        forAll(spField, fI)
        {
            bxyz.append(spField[fI]);
        }

        posLoc[myid]=bxyz;
        Pstream::allGatherList(posLoc);
        Field<GeoType> posGlob=ListListOps::combine<Field<GeoType>>(posLoc,accessOp<Field<GeoType>>());
        GeoTypeField["source"]=posGlob;
    }

    Field<GeoType> tpField(tP.size(),  pTraits<GeoType>::zero);

    if (ublist.size()>0)
    {
        Info<<"search boundary field..."<<endl;
        std::map<word,std::vector<label>>::iterator it;
        for (it=ublist.begin();it !=ublist.end();it++)
        {
            word pname=it->first;
            const std::vector<label> &blist=it->second;
            for (unsigned j=0;j<blist.size();j++)
            {
                label jc=blist[j];
                if (Pstream::parRun())
                {
                    tpField[j] = GeoTypeField["source"][jc];
                }
                else
                {
                    tpField[j] = spField[jc];
                }
            }
        }
    }

    // write abaqus inp file
    if (Pstream::master())
    {
        OFstream abaqusos
        (
            abaqusDir/
            "abaqusSurfaceData_" + word(field.name()) + "_" +
            field.mesh().time().timeName() + ".inp"
        );

        abaqusos << "**FOAM to ABAQUS load File" << endl;
        abaqusos << "**" << field.name() << endl;
        //abaqusos << "**"<< field.name() << " Units: "
        //         << field.dimensions() <<endl;
        abaqusos << "**Pressure Units: " << field.dimensions() <<endl;

        abaqusos<<"*DLOAD"<<endl;

        forAll(tP,faceID)
        {
            abaqusos
                << eid[faceID] << "," << field.name() << ","
                << tpField[faceID] << endl;
        }
    }

    sourceVTK.addFaceData(field.name(), spField);
    targetVTK.addFaceData(field.name(), tpField);
}

template<class GeoField>
void convertPressureToSIUnits
(
    GeoField& field
)
{
    // if incompressible, convert to [Pa]
    if (field.dimensions() == dimPressure/dimDensity)
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                field.time().constant(),
                field.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dimensionedScalar rho
        (
            transportProperties.lookup("rho")
        );

        field *= rho;
    }
}


template<class GeoField>
void print(const char* msg, Ostream& os, const UPtrList<GeoField>& flds)
{
    if (flds.size())
    {
        os  << msg;
        forAll(flds, i)
        {
            os  << ' ' << flds[i].name();
        }
        os  << endl;
    }
}



}

// Main program:
int main(int argc, char *argv[])
{
    //- user options
    argList::validArgs.clear();
    timeSelector::addOptions();
    argList::addBoolOption ("flipLoad","flip load sign");
    argList::addBoolOption ("flipNorm","flip Abaqus normals");

    argList::addBoolOption
    (
        "knnInterpolation", "knnInterpolation, default AMI"
    );

    argList::addOption
    (
        "translate",
        "vector",
        "translate by the specified <vector> - eg, '(1 0 0)'"
    );
    argList::addOption
    (
        "scale",
        "vector",
        "scale by the specified amount - eg, '(0.001 0.001 0.001)' for a "
        "uniform [mm] to [m] scaling"
    );
    argList::validOptions.insert
    (
        "targetSurface",
        "abaqus target surface in triSurface"
    );
    argList::validOptions.insert
    (
        "sourceSurface",
        "source patch from Foam or Elements case"
    );
    argList::validOptions.insert
    (
        "fields",
        "convert selected fields only. eg, '(p T U)' Only scalar vectors are "
        "supported"
    );

    argList args(argc, argv);

    Info<<"The foamToAbaqusSurfaceMap utility is running..."<<endl;

    //- target surface name
    word surfaceName;

    args.optionReadIfPresent("targetSurface", surfaceName);

    //- source patch, if specified
    word sourcePatchName;
    args.optionReadIfPresent("sourceSurface", sourcePatchName);
    bool isSourceSurfacesPresent (false);

    if (sourcePatchName.size()!=0)
    {
        isSourceSurfacesPresent = true;
        Info<<"The source is "<<sourcePatchName<<endl;
    }
    else
    {
        Info<<"The source is given by all the wall patches"<<endl;
    }

    wordHashSet selectedFields;
    const bool specifiedFields = args.optionReadIfPresent
    (
        "fields",
        selectedFields
    );

    bool knnInterpolation = args.optionFound("knnInterpolation");

    //flipDir
    bool invLoad = args.optionFound("flipLoad");

    scalar loadDir = 1;
    if (invLoad)
    {
        loadDir = -1;
        Info<<"pLoad direction flipped "<<loadDir<<endl;
    }

    //flipNorm
    bool flipAbaqusNorm = false;
    flipAbaqusNorm = args.optionFound("flipNorm");

#include "include/createTime.H"
#include "include/createMesh.H"

    //get time to work on different timesteps
    instantList timeDirs = timeSelector::select0(runTime, args);

    pointField targetPoints;
    faceList targetFaces;
    labelList eid;

    readINP(mesh, surfaceName, targetPoints, targetFaces, eid);

    vector v;

    if (args.optionReadIfPresent("translate", v))
    {
        Info<< "Translating points by " << v << endl;

        targetPoints += v;
    }

    if (args.optionReadIfPresent("scale", v))
    {
        Info<< "Scaling points by " << v << endl;

        targetPoints.replace(vector::X, v.x()*targetPoints.component(vector::X));
        targetPoints.replace(vector::Y, v.y()*targetPoints.component(vector::Y));
        targetPoints.replace(vector::Z, v.z()*targetPoints.component(vector::Z));
    }

    SubList<face> stargetFaces(targetFaces, targetFaces.size());
    primitivePatch targetPatch(stargetFaces, targetPoints);

    // calculate dimension
    label nWallFaces = 0;
    if (isSourceSurfacesPresent)
    {
        // use the specific patch as source
        label patchID(mesh.boundaryMesh().findPatchID(sourcePatchName));
        nWallFaces = mesh.boundary()[patchID].patch().size();
    }
    else
    {
        // group all the wall patches into the source
        forAll(mesh.boundaryMesh(),patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            nWallFaces = nWallFaces+mesh.boundary()[patchi].size();
        }
    }

    // fill wallFaces
    faceList wallFaces(nWallFaces);
    label faceCount = 0;
    if (isSourceSurfacesPresent)
    {
        // use the specific patch as source
        label patchID(mesh.boundaryMesh().findPatchID(sourcePatchName));
        forAll(mesh.boundary()[patchID],facei)
        {
            wallFaces[faceCount]= mesh.boundaryMesh()[patchID][facei];
            faceCount++;
        }
    }
    else
    {
        // group all the wall patches into the source
        forAll(mesh.boundaryMesh(),patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            forAll(mesh.boundary()[patchi],facei)
            {
                wallFaces[faceCount]= mesh.boundaryMesh()[patchi][facei];
                faceCount++;
            }
        }
    }

    // this is because primitivePatch doesn't like Lists
    SubList<face> sourceFaces(wallFaces, wallFaces.size());
    primitivePatch sourcePatch
    (
        sourceFaces,
        mesh.points()
    );

    // creating folder for results
    fileName abaqusDir = args.rootPath()/args.globalCaseName()/"AbaqusSurfaceData";
    if (Pstream::master())
    {
        if (isDir(abaqusDir))
        {
            rmDir(abaqusDir);
        }
        mkDir(abaqusDir);
    }

    autoPtr<AMIPatchToPatchInterpolation> amiInterpolatePtr(nullptr);
    std::map<word,KNN> bndMaps;
    std::map<word,std::vector<label>> faceIdMaps;


    if (!knnInterpolation)
    {
        //- AMI interpolation construction
        amiInterpolatePtr.reset
        (
            new AMIPatchToPatchInterpolation
            (
                sourcePatch,
                targetPatch,
                faceAreaIntersect::tmMesh,
                false,
                "partialFaceAreaWeightAMI",
                -1,
                flipAbaqusNorm //Abaqus normal control
            )
        );

    }
    else
    {
        constructKNNdata(sourcePatch, bndMaps, faceIdMaps);
    }

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobjectList objects(mesh, runTime.timeName());

        // Construct the vol fields
        // References the original mesh, but uses subsetted portion only.

        PtrList<volScalarField> vScalarFld;
        PtrList<volVectorField> vVectorFld;

        // Supported volume field types
        const wordHashSet vFieldTypes
        {
            volScalarField::typeName,
            volVectorField::typeName
        };

        if
        (
            candidateObjects
            (
                objects,
                vFieldTypes,
                specifiedFields,
                selectedFields
            ).size()
        )
        {
            readFields
            (
                mesh,
                objects,
                selectedFields,
                vScalarFld
            );
            print("    volScalar        :", Info, vScalarFld);

            readFields
            (
                mesh,
                objects,
                selectedFields,
                vVectorFld
            );
            print("    volVector        :", Info, vVectorFld);
        }

        simpleVTKWriter sourceVTK
        (
            sourceFaces,
            mesh.points()
        );
        sourceVTK.addFaceData("normalsField", sourcePatch.faceNormals());
        simpleVTKWriter targetVTK
        (
            targetFaces,
            targetPoints
        );
        targetVTK.addFaceData("normalsField", targetPatch.faceNormals());

        if (!knnInterpolation)
        {
            writeFieldsAMI
            (
                sourcePatch,
                targetPatch,
                isSourceSurfacesPresent,
                sourcePatchName,
                loadDir,
                amiInterpolatePtr,
                eid,
                vScalarFld,
                abaqusDir,
                sourceVTK,
                targetVTK
            );
            writeFieldsAMI
            (
                sourcePatch,
                targetPatch,
                isSourceSurfacesPresent,
                sourcePatchName,
                loadDir,
                amiInterpolatePtr,
                eid,
                vVectorFld,
                abaqusDir,
                sourceVTK,
                targetVTK
            );
        }
        else
        {
            writeFieldsKNN
            (
                sourcePatch,
                targetPatch,
                isSourceSurfacesPresent,
                sourcePatchName,
                loadDir,
                bndMaps,
                faceIdMaps,
                eid,
                vScalarFld,
                abaqusDir,
                sourceVTK,
                targetVTK
            );
            writeFieldsKNN
            (
                sourcePatch,
                targetPatch,
                isSourceSurfacesPresent,
                sourcePatchName,
                loadDir,
                bndMaps,
                faceIdMaps,
                eid,
                vVectorFld,
                abaqusDir,
                sourceVTK,
                targetVTK
            );
        }

        targetVTK.write
        (
            "target_"+mesh.time().timeName()+".vtk"
        );
        sourceVTK.write
        (
            "source_"+mesh.time().timeName()+".vtk"
        );
        //vtk print end
    }// end time loop 1

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
