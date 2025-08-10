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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2007 OpenCFD Ltd.
    (c) 2010-2012 ICON CFD

Application
    foamToNastran

Description
    Generate NASTRAN loads file.

\*---------------------------------------------------------------------------*/
#include "interpolations/primitivePatchInterpolation/PrimitivePatchInterpolation.H"
#include "cfdTools/general/include/fvCFD.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "primitives/random/Random/Random.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include <vector>
#if defined(WIN64) || defined(WIN32)
#include <math.h>
#endif
//#define DISPLAYERR 1
//define DISPLAYDIST 1

namespace Foam
{
    void writeVTKModified
    (
        const word& file_name,
        const pointField& allPoints,
        const DynamicList<face>& faces
    );
    void writeVTKModified
    (
        const word& file_name,
        const pointField& allPoints,
        const DynamicList<face>& faces,
        const scalarField& cellScalars
    );
    void writeVTKModified
    (
        const word& file_name,
        const pointField& allPoints,
        const DynamicList<face>& faces,
        const std::vector<scalar>& someScalars,
        const bool isWritePoints
    );
    void writeSTLModified
    (
        const word& file_name,
        const pointField& allPoints,
        const DynamicList<face>& faces
    );
    void readPrevPressureComputeError
    (
        const char* fileneame,
        std::vector<scalar>& cellScalars
    );
    void readPrevPointsComputeDist
    (
        const char* fileneame,
        const pointField& allPoints,
        std::vector<scalar>& nodalDist
    );
    bool lineTriangleIntersection
    (
        const Foam::triangle< Foam::point,
        const Foam::point& >& inputTri,
        const point& linept,
        const point& vect,
        point& ipoint
    );
    bool isSameClockDirection
    (
        const Foam::triangle< Foam::point, const Foam::point& >& inputTri,
        const Foam::point& edgeEnd
    );
 }


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Create patch from set of patches
autoPtr<indirectPrimitivePatch> makePatch
(
    const polyMesh& mesh,
    const labelList& patchIDs
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Count faces.
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        nFaces += pp.size();
    }

    // Collect faces.
    labelList addressing(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        label meshFaceI = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFaceI++;
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), addressing),
            mesh.points()
        )
    );
}


// Do weird things to extract number
static scalar parseNASCoord(const string& s)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign))());
        scalar exp = readScalar(IStringStream(s.substr(expSign+1))());

        if (s[expSign] == '-')
        {
            exp = -exp;
        }
        return mantissa*Foam::pow(scalar(10),exp);
    }
    else
    {
        return readScalar(IStringStream(s)());
    }
}


bool readNAS
(
    const IOobject& io,
    const fileName& OBJfileName,
    pointField& allPoints,
    DynamicList<face>& faces,
    DynamicList<label>& eid,
    DynamicList<label>& pid,
    const bool geoWrite=true
)
{
    autoPtr<IOobject> surfaceIO = io.clone();
    surfaceIO().rename(OBJfileName);

    IFstream OBJfile(typeFilePath<searchableSurface>(io));

    if (!OBJfile.good())
    {
        FatalErrorIn("readNAS(const fileName&)")
            << "Cannot read file " << surfaceIO().name()
            << exit(FatalError);
    }

    DynamicList<point> points;
    // Nastran index of point
    DynamicList<label> indices;

    // Done warnings per unrecognized command
    HashSet<word> unhandledCmd;

    while (OBJfile.good())
    {
        string line;
        OBJfile.getLine(line);

        if
        (line.size() == 0 || line[0] == '$' || line[0] == ' ' || line[0] == '+')
        {
            // Skip empty or comment
            continue;
        }

        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "CTRIA3")
        {
            face newFace;
            newFace.setSize(3);
            //label index, group, a, b, c;
            //lineStream >> index >> group >> a >> b >> c;
            label id = readLabel(IStringStream(line.substr(8,8))());
            label group = readLabel(IStringStream(line.substr(16,8))());
            newFace[0] = readLabel(IStringStream(line.substr(24,8))());
            newFace[1] = readLabel(IStringStream(line.substr(32,8))());
            newFace[2] = readLabel(IStringStream(line.substr(40,8))());

            faces.append(newFace);
            pid.append(group);
            eid.append(id);
        }
        else if (cmd == "CQUAD4")
        {
            face newFace;
            newFace.setSize(4);
            //label index, group, a, b, c, d;
            //lineStream >> index >> group >> a >> b >> c >> d;
            label id = readLabel(IStringStream(line.substr(8,8))());
            label group = readLabel(IStringStream(line.substr(16,8))());
            newFace[0] = readLabel(IStringStream(line.substr(24,8))());
            newFace[1] = readLabel(IStringStream(line.substr(32,8))());
            newFace[2] = readLabel(IStringStream(line.substr(40,8))());
            newFace[3] = readLabel(IStringStream(line.substr(48,8))());

            faces.append(newFace);
            pid.append(group);
            eid.append(id);
        }
        else if (cmd == "GRID")
        {
            //label index;
            //lineStream >> index;
            label index = readLabel(IStringStream(line.substr(8,8))());
            indices.append(index);

            scalar x = parseNASCoord(line.substr(24, 8));
            scalar y = parseNASCoord(line.substr(32, 8));
            scalar z = parseNASCoord(line.substr(40, 8));
            points.append(point(x, y, z));
        }
        else if (cmd == "GRID*")
        {
            // Assume on two lines with '*' continuation symbol on start of
            // second line. (comes out of Tgrid. Typical line (spaces truncated)
            // GRID*      126   0 -5.55999875E+02 -5.68730474E+02
            // *         2.14897901E+02
            string line2;
            OBJfile.getLine(line2);
            if (line2[0] != '*')
            {
                FatalErrorIn("readNAS(const fileName&)")
                    << "Expected continuation symbol '*' when reading GRID*"
                    << " (double precision coordinate) output by Tgrid" << nl
                    << "Read:" << line2 << nl
                    << "File:" << OBJfile.name()
                    << " line:" << OBJfile.lineNumber()
                    << exit(FatalError);
            }
            IStringStream lineStream(line.substr(10) + line2.substr(1));

            label index;
            lineStream >> index;
            indices.append(index);

            readScalar(lineStream); // What is this field?
            scalar x = readScalar(lineStream);
            scalar y = readScalar(lineStream);
            scalar z = readScalar(lineStream);
            points.append(point(x, y, z));
        }
        else if (unhandledCmd.insert(cmd))
        {
            // do nothing
        }
    }

    points.shrink();
    indices.shrink();
    faces.shrink();
    eid.shrink();
    pid.shrink();

    {
        // Build inverse mapping (index to point)
        Map<label> indexToPoint(2*indices.size());
        forAll(indices, i)
        {
            indexToPoint.insert(indices[i], i);
        }

        // Relabel triangles
        forAll(faces, i)
        {
            face& f = faces[i];
            forAll(f,fp)
            {
                f[fp] = indexToPoint[f[fp]];
            }
        }
    }

    // Transfer DynamicLists to straight ones.
    allPoints.transfer(points);
    points.clear();

    if (geoWrite) {
       writeVTKModified ( surfaceIO().name() , allPoints, faces );
       //writeSTLModified ( surfaceIO().name() , allPoints, faces );
   }

    return true;
}

namespace Foam
{
void writeVTKModified
(
    const word& file_name,
    const pointField& allPoints,
    const DynamicList<face>& faces
)
{
    // 0. opening data output file in ascii format
      char outfile [512U];

      strcpy( outfile, file_name.c_str() );
      strcat( outfile, ".vtk" );

      ofstream ofs;
      ofs.open( outfile );
      if (!ofs)
      {
           Info<<" foamToNastran::writeVTKModified(): Output file "
                <<outfile<<" could not be opened."<< endl;
           return;
      }

       ofs.precision(10U);

       // 1. writing the file header
       // --------------------------
       ofs <<"# vtk DataFile Version 2.0"<< std::endl;
       ofs <<"Finite-element dataset: OpenFOAM variable: dummy"<< std::endl;
       ofs <<"ASCII"<< std::endl << std::endl;

       const label nnodes( allPoints.size() );

       // 2. writing node coordinates
       // ---------------------------
       ofs <<"DATASET UNSTRUCTURED_GRID"<< std::endl;
       ofs <<"POINTS " << nnodes <<" float"<< std::endl;

       forAll(allPoints, i )
       {
           ofs<<allPoints[i][0]<<" "<<allPoints[i][1]<<" "<<allPoints[i][2]
              <<std::endl;
       }
       ofs << std::endl;

       // 3. writing CELLS (cell-size and member nodes (point)) -corect
       // Cells
       // -----------------------------------------------------
       label total_entries(0);
       forAll( faces, i )
       {
         total_entries+=( faces[i].size()+1 );
       }

       const label nelements ( faces.size() );

      ofs <<"CELLS "<< nelements <<" "<< total_entries << std::endl;
      forAll( faces, i )
      {
         ofs << faces[i].size();
         forAll( faces[i],j) ofs<<" "<<faces[i][j];
         ofs << std::endl;
       }
       ofs << std::endl;
       // 4. writing CELL_TYPES - for the
       // ---------------------
       ofs <<"CELL_TYPES "<< nelements << std::endl;
       forAll( faces, i )
       {
           if (faces[i].size()==3U)
           {
               ofs<<"5";
               if ((i%10) || !i) ofs<<" ";
               else ofs<<std::endl;
           }

           if (faces[i].size()==4U)
           {
               ofs<<"9";
               if ((i%10) || !i) ofs<<" ";
               else ofs<<std::endl;

           }
       }
       ofs << std::endl;

       ofs.close();
       Info<<" foamToNastran::writeVTKModified(): VTK file"<<outfile
           <<" sucessfuly written."<<endl;

} // writeVTKModified

//
//Overloaded with one scalar field output
void writeVTKModified
(
    const word& file_name,
    const pointField& allPoints,
    const DynamicList<face>& faces,
    const scalarField& cellScalars
)
{
    // 0. opening data output file in ascii format
      char outfile [512U];

      strcpy( outfile, file_name.c_str() );
      strcat( outfile, ".vtk" );

      ofstream ofs;
      ofs.open( outfile );
      if (!ofs)
      {
           Info<<" foamToNastran::writeVTKModified(): Output file "
                <<outfile<<" could not be opened."<< endl;
           return;
      }

       ofs.precision(10U);

       // 1. writing the file header
       // --------------------------
       ofs <<"# vtk DataFile Version 2.0"<< std::endl;
       ofs <<"Finite-element dataset: OpenFOAM variable: dummy"<< std::endl;
       ofs <<"ASCII"<< std::endl << std::endl;

       const label nnodes( allPoints.size() );

       // 2. writing node coordinates
       // ---------------------------
       ofs <<"DATASET UNSTRUCTURED_GRID"<< std::endl;
       ofs <<"POINTS " << nnodes <<" float"<< std::endl;

       forAll(allPoints, i )
       {
           ofs<<allPoints[i][0]<<" "<<allPoints[i][1]<<" "
              <<allPoints[i][2]<<std::endl;
       }
       ofs << std::endl;

       // 3. writing CELLS (cell-size and member nodes (point)) -corect
       // Cells
       // -----------------------------------------------------
       label total_entries(0);
       forAll( faces, i )
       {
           total_entries+=( faces[i].size()+1 );
       }

       const label nelements ( faces.size() );

       ofs <<"CELLS "<< nelements <<" "<< total_entries << std::endl;
       forAll( faces, i )
       {
           ofs << faces[i].size();
           forAll( faces[i],j) ofs<<" "<<faces[i][j];
           ofs << std::endl;
       }
       ofs << std::endl;
       // 4. writing CELL_TYPES - for the
       // ---------------------
       ofs <<"CELL_TYPES "<< nelements << std::endl;
       forAll( faces, i )
       {
           if (faces[i].size()==3U)
           {
               ofs<<"5";
               if ((i%10) || !i) ofs<<" ";
               else ofs<<std::endl;
           }

           if (faces[i].size()==4U)
           {
               ofs<<"9";
               if ((i%10) || !i) ofs<<" ";
               else ofs<<std::endl;
           }
       }
       ofs<<std::endl;

       if (cellScalars.size())
       {
           if (cellScalars.size() !=nelements)
           {
               std::cerr<<" ***ERROR:foamToNastran::writeVTKModified()(2): "
                        <<"wrong dimension of provided scalar field ="
                        <<cellScalars.size()<<" vs ele="<<nelements<<std::endl;
               return;
           }
           ofs<<"CELL_DATA "<<nelements<<std::endl;
           ofs<<"FIELD attributes 1"<<std::endl;
           ofs<<"value 1 "<<nelements<<" float"<<std::endl;
           forAll( cellScalars, i )
           {
               ofs<<cellScalars[i];
               if ((i%10) || !i) ofs<<" ";
               else ofs<<std::endl;
           }
       }

       ofs << std::endl;
       ofs.close();
       Info<<" foamToNastran::writeVTKModified(): VTK file"<<outfile
           <<" sucessfuly written."<<endl;
} // writeVTKModified overloaded FOAM

//Std Overloaded with one scalar field output
void writeVTKModified
(
    const word& file_name,
    const pointField& allPoints,
    const DynamicList<face>& faces,
    const std::vector<scalar>& someScalars,
    const bool isWritePoints
)
{
    // 0. opening data output file in ascii format
      char outfile [512U];

      strcpy( outfile, file_name.c_str() );
      strcat( outfile, ".vtk" );

      ofstream ofs;
      ofs.open( outfile );
      if (!ofs)
      {
           Info<<" foamToNastran::writeVTKModified(): Output file "<<outfile
                <<" could not be opened."<< endl;
           return;
      }

      ofs.precision(10U);

       // 1. writing the file header
       // --------------------------
       ofs <<"# vtk DataFile Version 2.0"<< std::endl;
       ofs <<"Finite-element dataset: OpenFOAM variable: dummy"<< std::endl;
       ofs <<"ASCII"<< std::endl << std::endl;

       const label nnodes( allPoints.size() );

       // 2. writing node coordinates
       // ---------------------------
       ofs <<"DATASET UNSTRUCTURED_GRID"<< std::endl;
       ofs <<"POINTS " << nnodes <<" float"<< std::endl;

       forAll(allPoints, i )
       {
           ofs<<allPoints[i][0]<<" "<<allPoints[i][1]<<" "<<allPoints[i][2]
              <<std::endl;
       }
       ofs << std::endl;

       // 3. writing CELLS (cell-size and member nodes (point)) -corect
       // Cells
       // -----------------------------------------------------
       label total_entries(0);
       forAll( faces, i )
       {
           total_entries+=( faces[i].size()+1 );
       }

       const unsigned int nelements ( faces.size() );

      ofs <<"CELLS "<< nelements <<" "<< total_entries << std::endl;
      forAll( faces, i )
      {
          ofs << faces[i].size();
          forAll( faces[i],j) ofs<<" "<<faces[i][j];
          ofs << std::endl;
      }
      ofs << std::endl;
      // 4. writing CELL_TYPES - for the
      // ---------------------
      ofs <<"CELL_TYPES "<< nelements << std::endl;
      forAll( faces, i )
      {
          if (faces[i].size()==3U)
          {
              ofs<<"5";
              if ((i%10)) ofs<<" ";
              else ofs<<std::endl;
          }

          if (faces[i].size()==4U)
          {
              ofs<<"9";
              if ((i%10) || !i) ofs<<" ";
              else ofs<<std::endl;
          }
      }// faces

      ofs << std::endl;

      if (isWritePoints)
      {
          if (someScalars.size())
          {
              if (someScalars.size() !=static_cast <unsigned int> (nnodes))
              {
                  std::cerr<<" ***ERROR:foamToNastran::writeVTKModified()(3): "
                           <<"wrong dimension of provided scalar field ="
                           <<someScalars.size()<<" vs ele="<<nnodes<<std::endl;
                  return;
              }
              ofs<<std::endl<<"POINT_DATA "<<nnodes<<std::endl;
              ofs<<"FIELD attributes 1"<<std::endl;
              ofs<<"value 1 "<<nnodes<<" float"<<std::endl;
              for (unsigned int i=0; i<someScalars.size(); i++)
              {
                  ofs<<someScalars[i];
                  if ((i%10) || !i) ofs<<" ";
                  else ofs<<std::endl;
              }
          } //if size
      }
      else
      {
          if (someScalars.size())
          {
              if (someScalars.size() !=nelements)
              {
                  std::cerr<<" ***ERROR:foamToNastran::writeVTKModified()(3): "
                           <<"wrong dimension of provided scalar field ="
                           <<someScalars.size()<<" vs ele="<<nelements
                           <<std::endl;
                  return;
              }
              ofs<<std::endl<<"CELL_DATA "<<nelements<<std::endl;
              ofs<<"FIELD attributes 1"<<std::endl;
              ofs<<"value 1 "<<nelements<<" float"<<std::endl;
              for (unsigned int i=0; i<someScalars.size(); i++)
              {
                  ofs<<someScalars[i];
                  if ((i%10) || !i) ofs<<" ";
                  else ofs<<std::endl;
              }
          } //if size
      } //else

       ofs << std::endl;
       ofs.close();
       Info<<" foamToNastran::writeVTKModified(): VTK file"<<outfile
           <<" sucessfuly written."<<endl;

} // writeVTKModified overloaded


void writeSTLModified
(
    const word& file_name,
    const pointField& allPoints,
    const DynamicList<face>& faces
)
{
    // 0. opening data output file in ascii format
      char outfile [512U];

      strcpy( outfile, file_name.c_str() );
      strcat( outfile, ".stl" );

      ofstream ofs;
      ofs.open( outfile );
      if (!ofs)
      {
           Info<<" foamToNastran::writeSTLModified(): Output file "<<outfile
                <<" could not be opened."<< endl;
           return;
      }

      ofs.precision(10U);

      ofs<<"solid nastran_target_surface"<<std::endl;

      const label nelements ( faces.size() );
      Info<<" foamToNastran::writeSTLModified(): saving "<<nelements
           <<" surface elements to file "<<outfile<<endl;

      forAll( faces, i )
      {
          if (faces[i].size()==3)
          {
              Foam::point a(allPoints[ faces[i][0] ]);
              Foam::point b(allPoints[ faces[i][1] ]);
              Foam::point c(allPoints[ faces[i][2] ]);
              Foam::triangle< point, const point& > tri(a, b, c);

              Foam::vector vec(tri.areaNormal());

              ofs<<" facet normal "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]
                 <<std::endl;
              ofs<<"  outer loop"<<std::endl;
              ofs<<"   vertex "<<a[0]<<" "<<a[1]<<" "<<a[2]<<std::endl;
              ofs<<"   vertex "<<b[0]<<" "<<b[1]<<" "<<b[2]<<std::endl;
              ofs<<"   vertex "<<c[0]<<" "<<c[1]<<" "<<c[2]<<std::endl;
              ofs<<"  endloop"<<std::endl;
              ofs<<" endfacet"<<std::endl;
          }
          else if (faces[i].size()==4)
          {
              Foam::point a(allPoints[ faces[i][0] ] );
              Foam::point b(allPoints[ faces[i][1] ] );
              Foam::point c(allPoints[ faces[i][2] ] );
              Foam::point d(allPoints[ faces[i][3] ]);
              Foam::triangle< point, const point& > tri1(a, b, c);
              Foam::triangle< point, const point& > tri2(a, c, d);

              Foam::vector vec1(tri1.areaNormal());
              Foam::vector vec2(tri2.areaNormal());

            ofs<<" facet normal "<<vec1[0]<<" "<<vec1[1]<<" "<<vec1[2]
               <<std::endl;
            ofs<<"  outer loop"<<std::endl;
            ofs<<"   vertex "<<a[0]<<" "<<a[1]<<" "<<a[2]<<std::endl;
            ofs<<"   vertex "<<b[0]<<" "<<b[1]<<" "<<b[2]<<std::endl;
            ofs<<"   vertex "<<c[0]<<" "<<c[1]<<" "<<c[2]<<std::endl;
            ofs<<"  endloop"<<std::endl;
            ofs<<" endfacet"<<std::endl;

            ofs<<" facet normal "<<vec2[0]<<" "<<vec2[1]<<" "<<vec2[2]
               <<std::endl;
            ofs<<"  outer loop"<<std::endl;
            ofs<<"   vertex "<<a[0]<<" "<<a[1]<<" "<<a[2]<<std::endl;
            ofs<<"   vertex "<<c[0]<<" "<<c[1]<<" "<<c[2]<<std::endl;
            ofs<<"   vertex "<<d[0]<<" "<<d[1]<<" "<<d[2]<<std::endl;
            ofs<<"  endloop"<<std::endl;
            ofs<<" endfacet"<<std::endl;

          }
          else
              std::cerr<<" ***ERROR:foamToNastran::writeSTLModified(): "
                       <<"unsupported surface element with "<<faces[i].size()
                       <<" nodes!"
                       <<std::endl;
      }
       ofs<<"endsolid nastran_target_surface"<<std::endl;
       ofs << std::endl;

       ofs.close();
       Info<<" foamToNastran::writeSTLModified(): STL file"<<outfile
           <<" sucessfuly written."<<endl;

} // writeSTLModified


void readPrevPressureComputeError
(
    const char* fileneame,
    std::vector<scalar>& cellScalars
)
{
    // 0. opening data output file in ascii format
      char infile [512U];
      strcpy( infile, fileneame );
      strcat( infile, ".dat" );

      ifstream ifs;
      ifs.open( infile );
      if (!ifs)
      {
           std::cerr <<" foamToNastran::readPrevPressureComputeError(): Input "
                     <<"file "<<infile<<" could not be opened."<< std::endl;
           return;
      }

      if (!cellScalars.empty())
          cellScalars.erase( cellScalars.begin(),cellScalars.end() );


      unsigned int nelements(0);
      ifs>>nelements;

      for (unsigned int i=0; i<nelements; i++)
      {
          scalar curr_value;
          ifs>>curr_value;
          cellScalars.push_back(curr_value);
      }

      ifs.close();
      return;
}


scalar computeHermitianSpectra
(
    const pointField& allPoints,
    const DynamicList<face>& faces
)
{
   scalar dSpectra (1.0);

   /*
   //scalar sSpectra (0.0);

   scalar affineRatio(0.0);
   forAll( faces, i )
   {
       vector  aVect(1.0,1.0,1.0);
       scalarField v(faces[i].size(),0.866);
       forAll( faces[i], j )
       {
           aVect=allPoints[ faces[i][j] ];
       }
       affineRatio+=( faces[i].areaInContact( allPoints,v ) );
   }

   affineRatio/=faces.size();

   //sSpectra =affineRatio;
   */

   return dSpectra;
}

void readPrevPointsComputeDist
(
    const char* fileneame,
    const pointField& allPoints,
    std::vector<scalar>& nodalSclars
)
{
      // 0. opening data output file in ascii format
      char infile [512U];
      strcpy( infile, fileneame );
      strcat( infile, ".dat" );

      ifstream ifs;
      ifs.open( infile );
      if (!ifs)
      {
          std::cerr <<" foamToNastran::readPrevPointsComputeDist(): "<<
              "Input file "<<infile<<" could not be opened."<< std::endl;
          return;
      }

      if (!nodalSclars.empty())
          nodalSclars.erase( nodalSclars.begin(),nodalSclars.end() );


      unsigned int nnodes(0);
      ifs>>nnodes;

      if (nnodes!= static_cast <unsigned int> (allPoints.size()))
      {
          std::cerr <<" foamToNastran::readPrevPointsComputeDist(): point "
                    <<"arrays do NOT match: in file="<<nnodes<<", in field="
                    <<allPoints.size()<< std::endl;
          return;
      }

      for (unsigned int i=0; i<nnodes; i++)
      {
          scalar x,y,z;
          ifs>>x;
          ifs>>y;
          ifs>>z;
          scalar dist(0.0);

          dist=sqrt( (x-allPoints[i][0])*(x-allPoints[i][0]) +
                     (y-allPoints[i][1])*(y-allPoints[i][1]) +
                     (z-allPoints[i][2])*(z-allPoints[i][2]) );
          nodalSclars.push_back(dist);
      }

      ifs.close();
      return;
}


bool lineTriangleIntersection
(
    const Foam::triangle< Foam::point, const Foam::point& >& inputTri,
    const Foam::point& linept,
    const Foam::point& vect,
    Foam::point& ipoint
)
{
  const vector triNormal(inputTri.areaNormal());

  const scalar dotProduct
      (triNormal[0]*vect[0]+triNormal[1]*vect[1]+triNormal[2]*vect[2]);

  scalar tolerance(1e-12);

  if (dotProduct<tolerance)
  {
      const scalar time ( -(triNormal[0]*(linept[0]-inputTri.a()[0])
                            +triNormal[1]*(linept[1]-inputTri.a()[1])+
                            triNormal[2]*(linept[2]-inputTri.a()[2]))/
                          (triNormal[0]*vect[0]+triNormal[1]*vect[1]
                           +triNormal[2]*vect[2]) );

      // Here ther is a lot of problems, related what is the orientation,
      // what starts where, etc.
      //if(time<0.0) return false;

      ipoint[0]=linept[0] + vect[0]*time;
      ipoint[1]=linept[1] + vect[1]*time;
      ipoint[2]=linept[2] + vect[2]*time;

      Foam::triangle< point, const point& >
          prom0(inputTri.a(),inputTri.b(),ipoint);
      Foam::triangle< point, const point& >
          prom1(inputTri.b(),inputTri.c(),ipoint);
      Foam::triangle< point, const point& >
          prom2(inputTri.c(),inputTri.a(),ipoint);

      if (isSameClockDirection(prom0,triNormal))
      {
          if (isSameClockDirection(prom1,triNormal))
          {
              if (isSameClockDirection(prom2,triNormal))
              {
                  return true;
              }
          }
      }
  }
  return false;
}



bool  isSameClockDirection
(
    const Foam::triangle< point, const point& >& inputTri,
    const point& edgeEnd
)
{
    const vector triNormal(inputTri.areaNormal());

    if
    (
        (triNormal[0]*edgeEnd[0]
         + triNormal[1]*edgeEnd[1]
         + triNormal[2]*edgeEnd[2] )<0.0
     )
        return false;
    else
        return true;
}



} // namespace Foam


// Main program:

int main(int argc, char *argv[])
{

    const bool ldebug (false);

//    argList::validArgs.append("output format");
#include "include/addTimeOptions.H"
#include "include/setRootCase.H"
#include "include/createTime.H"
#include "include/createMesh.H"

    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#include "include/checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);

    // Read NASTRANS dictionary
    IOdictionary nastranDict
    (
       IOobject
       (
            "nastranDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );

    // whether using unsteady data
    bool unsteady = true;
    if (nastranDict.found("unsteady"))
    {
        unsteady =  readBool(nastranDict.lookup("unsteady"));
    }
    else
    {
        Info<<"foamToNastran::main(): No unsteady keyword in nastranDict. "
            <<"Using Pmean"<<endl;
    }
    word depVar;
    if (unsteady)
    {
        depVar = "pMean";
    }
    else
    {
        depVar = "p";
    }

    const bool geoout = nastranDict.lookupOrDefault<Switch>("geooutput",false);

    // Scaling factor for unit handling
    const scalar pressureUnitsScaling =
        nastranDict.lookupOrDefault<scalar>("pressureUnitsScaling",1.0);

    const scalar pressurePhysicsScaling =
        nastranDict.lookupOrDefault<scalar>("pressurePhysicsScaling",1.0);

    label maxIter =
        nastranDict.lookupOrDefault<label>("maxProjectionIterations",5);

    scalar scalingConstant =
        nastranDict.lookupOrDefault<scalar>("scalingConstant",2.0);


    // Find surfaces TO which we want to perform mapping surface==TARGET
    PtrList<dictionary> targetDicts = nastranDict.lookup("targets");
    fileNameList    targetFileNames(targetDicts.size());
    List<labelList> targetOutputPIDs(targetDicts.size());
    List<scalar>    targetBackPressure(targetDicts.size());

    scalar globalBackPressure =
        nastranDict.lookupOrDefault<scalar>("globalBackPressure", 0.0);

    forAll(targetDicts, i)
    {
        const dictionary& dict = targetDicts[i];
        targetFileNames[i] = fileName(dict.lookup("file"));
        targetOutputPIDs[i] = labelList(dict.lookup("pids"));
        if (dict.found("backPressure"))
        {
           targetBackPressure[i] = readScalar(dict.lookup("backPressure"));
        }
        else
        {
           targetBackPressure[i] = globalBackPressure;
        }
    }


    // Find surfaces FROM which we want to perform mapping
    PtrList<dictionary> sourceDicts = nastranDict.lookup("sources");
    fileNameList sourceFileNames( sourceDicts.size() );

    bool isSourceSurfacesPresent (false);
    if (sourceDicts.size()) isSourceSurfacesPresent=true;

    forAll(sourceDicts, i)
    {
        const dictionary& dict = sourceDicts[i];
        sourceFileNames[i] = fileName(dict.lookup("file"));
    }

    if (isSourceSurfacesPresent)
    {
       if (sourceDicts.size() != targetDicts.size())
           FatalError << "Inconsistent number of NASTRAN source and target "
                      <<"projection files"<< exit(FatalError);
    }

    List<pointField> nastranTargetFaceCentres(targetFileNames.size());
    List<labelList> nastranTargetEID(targetFileNames.size());
    List<labelList> nastranTargetPID(targetFileNames.size());

    // All points - global accounting
    pointField allPointsGlobal;
    // Faces in terms of Nastran point indices
    DynamicList<face> facesGlobal;

    //All related to TARGET
    pointField allPointsTargetGlobal;
    // Faces in terms of Nastran point indices
    DynamicList<face> facesTargetGlobal;

    // Load all surfaces specified in autoHexMeshDict.
    forAll(targetFileNames, i)
    {
        word ext = targetFileNames[i].ext();
        if (ext == "nas")
        {
            // Nastran index of elements
            DynamicList<label> eid;
            // Nastran index of patches
            DynamicList<label> pid;
            // coordinates of point
            pointField allPoints;
            // Faces in terms of Nastran point indices
            DynamicList<face> faces;

            readNAS
            (
                IOobject
                (
                    "",                                 // dummy name
                    runTime.constant(),                 // directory
                    "triSurface",                       // instance
                    runTime,                            // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                targetFileNames[i],
                allPoints,
                faces,
                eid,
                pid,
                geoout
            );


            allPointsTargetGlobal=allPoints;
            facesTargetGlobal=faces;

            pointField fc(faces.size());
            forAll(fc, faceI)
            {
                fc[faceI] = faces[faceI].centre(allPoints);
            }
            nastranTargetFaceCentres[i] = fc;
            nastranTargetEID[i] = labelList(eid);
            nastranTargetPID[i] = labelList(pid);
        }
        else
        {
            FatalError << "File '"
                       << targetFileNames[i] << "' is not a "
                       << "nastran format file."
                       << exit(FatalError);
        }
    }


    //Process nastaran SOURCE surface  (if present)
    List<pointField> nastranSourceFaceCentres(targetFileNames.size());
    List<labelList> nastranSourceEID(targetFileNames.size());
    List<labelList> nastranSourcePID(targetFileNames.size());

    // All points
    pointField allPointsSourceGlobal;
    // Faces in terms of Nastran point indices
    DynamicList<face> facesSourceGlobal;


    // Load all surfaces specified in autoHexMeshDict.
    forAll(sourceFileNames, i)
    {
        word ext = sourceFileNames[i].ext();
        if (ext == "nas")
        {
            // Nastran index of elements
            DynamicList<label> eidl;
            // Nastran index of patches
            DynamicList<label> pidl;
            // coordinates of point
            pointField allPointsLocal;
            // Faces in terms of Nastran point indices
            DynamicList<face> facesl;

            readNAS
            (
                IOobject
                (
                    "",                                 // dummy name
                    runTime.constant(),                 // directory
                    "triSurface",                       // instance
                    runTime,                            // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                sourceFileNames[i],
                allPointsLocal,
                facesl,
                eidl,
                pidl,
                geoout
            );

            allPointsSourceGlobal=allPointsLocal;
            facesSourceGlobal=facesl;

            pointField fc(facesl.size());
            forAll(fc, faceI)
            {
                fc[faceI] = facesl[faceI].centre(allPointsLocal);
            }
            nastranSourceFaceCentres[i] = fc;
            nastranSourceEID[i] = labelList(eidl);
            nastranSourcePID[i] = labelList(pidl);
        }
        else
        {
            FatalError << "Sorce surface file '"
                       << targetFileNames[i] << "' is not a "
                       << "nastran format file."
                       << exit(FatalError);
        }
    }


    // Construct table of patches to include.
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    DynamicList<label> includePatches(bMesh.size());
    for (label patchI = 0; patchI < bMesh.size(); patchI++)
    {
        if (!isA<processorPolyPatch>(bMesh[patchI]))
        {
            includePatches.append(patchI);
        }
    }
    includePatches.shrink();

    // Build single patch out of included patches
    indirectPrimitivePatch pp = makePatch(mesh, includePatches);

   // Setup indirectPrimitivePatch interpolation
    PrimitivePatchInterpolation<indirectPrimitivePatch> ppInterp(pp);

    // Create triSurface from boundary mesh
    triSurface localSurface
    (
        triSurfaceTools::triangulate
        (
            mesh.boundaryMesh(),
            labelHashSet(includePatches)
        )
    );

    // Random number generator
    Random rndGen(354543);

    // bb of surface
    treeBoundBox bb(localSurface.localPoints());

    // Calculate an approximate length scale for the mesh
    scalar maxEdgeLength = bb.mag();
    maxEdgeLength = returnReduce(maxEdgeLength, sumOp<scalar>());

    scalar maxEdgeLengthSqr = maxEdgeLength * maxEdgeLength;

    triSurfaceSearch searchSelectSurf
    (
        localSurface,
        indexedOctree<treeDataTriSurface>::perturbTol(),
        8
     );

    // search engine
    const indexedOctree<treeDataTriSurface>& surfaceTree =
        searchSelectSurf.tree();

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho
    (
        "rho", dimDensity, transportProperties.lookup("rho")
    );


    //Choosing where to project
    List<pointField> nastranFaceCentres;
    List<labelList> nastranEID;
    List<labelList> nastranPID;

    // All points
    //pointField allPoints;
    // Faces in terms of Nastran point indices
    //DynamicList<face> allFaces;


    PtrList<dictionary> surfaceDicts = nastranDict.lookup("targets");
    fileNameList    surfaceFileNames   (surfaceDicts.size());
    List<labelList> surfaceOutputPIDs  (surfaceDicts.size());
    List<scalar>    surfaceBackPressure(surfaceDicts.size());


    if (isSourceSurfacesPresent)
    {
        Info<<endl<<" WARNING: foamToNastran::main(): perfoming SURFACE "
            <<"to SURFACE projection..."<<endl;

        nastranFaceCentres=nastranSourceFaceCentres;
        nastranEID=nastranSourceEID;
        nastranPID=nastranSourcePID;
        allPointsGlobal=allPointsSourceGlobal;
        facesGlobal=facesSourceGlobal;
        surfaceDicts=sourceDicts;
        surfaceFileNames=targetFileNames;
    }
    else
    {
        nastranFaceCentres=nastranTargetFaceCentres;
        nastranEID=nastranTargetEID;
        nastranPID=nastranTargetPID;
        allPointsGlobal=allPointsTargetGlobal;
        facesGlobal=facesTargetGlobal;
        surfaceDicts=targetDicts;
        surfaceFileNames=sourceFileNames;
    }

    forAll(Times, timeI)
    {
        runTime.setTime(Times[timeI], timeI);

        Info<< "Time " << Times[timeI].name() << endl;

        volScalarField pMean
        (
            IOobject
            (
                depVar,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        if (pMean.dimensions() != dimPressure)
        {
            pMean *= rho;
        }

        scalarField patchPMean(pp.size());
        forAll(patchPMean, faceI)
        {
            label patchID =
                mesh.boundaryMesh().whichPatch(pp.addressing()[faceI]);
            label pFaceI =
                pp.addressing()[faceI] - mesh.boundaryMesh()[patchID].start();
            patchPMean[faceI] = pMean.boundaryField()[patchID][pFaceI];
        }

        // Interpolate from face centres to points
        scalarField patchPointsPMean( ppInterp.faceToPointInterpolate(patchPMean) );

        label forceNumber = 1;
        forceNumber +=
            nastranDict.lookupOrDefault<label>("forceNumberOffset",0);

        label nFailedAttempts(0);
        scalar hermitianSpectra (0.866);

        forAll(surfaceFileNames, surfI)
        {
            const pointField& nastranCentres = nastranFaceCentres[surfI];
            scalarField pInterpolated(nastranCentres.size(), 0.0);
            scalarField closestDistance(nastranCentres.size(), GREAT);
            forAll(nastranCentres, faceI)
            {
                const point fc(nastranCentres[faceI]);

                bool projectionFound (false);

                pointIndexHit hit
                (
                    surfaceTree.findNearest(fc, maxEdgeLengthSqr)
                );

                if (hit.hit()) projectionFound=true;

                label niter(0);

                while (!projectionFound && niter<maxIter)
                {
                    maxEdgeLengthSqr*=scalingConstant;
                    hit=  surfaceTree.findNearest(fc, maxEdgeLengthSqr);
                    if (hit.hit()) { projectionFound=true; }
                    niter++;
                    if (ldebug)
                        Pout<<"foamToNastran::main(): foamToNastran(): "
                            <<"projection iteration="<<niter<<endl;
                }

                scalar distHit (GREAT);

                if (projectionFound)
                {
                    distHit = mag(hit.hitPoint() - fc);

                   // Perform inverse distance weighted interpolation
                   label index = hit.index();
                   const labelledTri& f = localSurface.localFaces()[index];
                   scalarField weights(f.size());

                   forAll(f, fI)
                    {
                        const point& tPoint =  localSurface.points()[f[fI]];
                        scalar dist = mag(tPoint - fc);
                        weights[fI] = 1.0 / (dist + VSMALL);
                        pInterpolated[faceI] +=
                            weights[fI] * patchPointsPMean[f[fI]];
                    }
                   pInterpolated[faceI] /= sum(weights);
                   closestDistance[faceI] = distHit;
                }
                else
                {
                    if (localSurface.points().size() != 0)
                    {
                        // No closest point found, so set to mapped pressure
                        //to zero
                        Pout<<"foamToNastran::main(): No closest mesh surface "
                            <<"found to NASTRAN face centre: " << fc
                            <<" on surface: "<< targetFileNames[surfI].lessExt()
                            <<" .Setting mapped pressure to zero"<<endl;
                    }
                    pInterpolated[faceI] = 0.0;
                    nFailedAttempts++;
                }

            }
            //calculate minimum closest distance over all processors
            scalarField minClosestDistance = closestDistance;

            Pstream::listCombineReduce(minClosestDistance, minOp<scalar>());

            labelList foundClosest(nastranCentres.size(), 0);
            forAll(minClosestDistance, faceI)
            {
                if (closestDistance[faceI] > 1.01 * minClosestDistance[faceI])
                {
                    // set mapped pressure to zero if not closest
                    pInterpolated[faceI] = 0.0;
                }
                else
                {
                    foundClosest[faceI]++;
                }
            }

            Pstream::listCombineReduce(pInterpolated, plusOp<scalar>());
            Pstream::listCombineReduce(foundClosest, plusOp<label>());

            forAll(pInterpolated, faceI)
            {
                if (foundClosest[faceI] > 0)
                {
                    pInterpolated[faceI] /= foundClosest[faceI];
                }
            }

            // Output ploads file

            if (Pstream::master() || !Pstream::parRun())
            {
                std::vector<scalar> cellValuesVector;
                forAll(targetOutputPIDs[surfI], i)
                {
                    label pid = targetOutputPIDs[surfI][i];
                    fileName ploadFileName = targetFileNames[surfI].lessExt()
                                           +"_pid"+name(pid)
                                           +"_time"+Times[timeI].name()
                                           +"."+"pload";
                    if (ldebug) Info<<"Writing PLOAD file: "<<ploadFileName<<endl;
                    if (Pstream::parRun())
                    {
                        ploadFileName
                            =  runTime.rootPath()/runTime.caseName()/".."
                            /"PLOAD"/ploadFileName;
                    }
                    else
                    {
                        ploadFileName
                            = runTime.rootPath()/runTime.caseName()/"PLOAD"
                            /ploadFileName;
                    }

                    if (Pstream::parRun())
                    {
                        if (!Foam::isDir(runTime.rootPath()/runTime.caseName()
                                      /".."/"PLOAD"))
                        {
                            Foam::mkDir(runTime.rootPath()/runTime.caseName()
                                        /".."/"PLOAD");
                        }
                    }
                    else
                    {
                        if (!Foam::isDir(runTime.rootPath()/runTime.caseName()
                                      /"PLOAD"))
                        {
                            Foam::mkDir(runTime.rootPath()/runTime.caseName()
                                        /"PLOAD");
                        }
                    }
                    autoPtr<OFstream> ploadFilePtr(new OFstream(ploadFileName));
                    label numberWithPID (0);
                    forAll(pInterpolated, faceI)
                    {
                        if (nastranPID[surfI][faceI] == pid)
                        {
                            scalar p = pressureUnitsScaling
                                * pressurePhysicsScaling
                                * (pInterpolated[faceI] +
                                   targetBackPressure[surfI]);

                            ploadFilePtr().unsetf(ios_base::fixed
                                                  | ios_base::scientific);


                            if (p > 0)
                            {
                                ploadFilePtr().precision(3);
                            }
                            else
                            {
                                ploadFilePtr().precision(2);
                            }

                            ploadFilePtr()<< setw(8) <<"PLOAD2  "
                                          << setw(7) <<forceNumber
                                          << setw(1) <<" "
                                          << setw(8) <<p
                                          << setf(ios_base::left)
                                          << setw(8) <<nastranEID[surfI][faceI]
                                          << endl;
                            ploadFilePtr().unsetf(ios_base::left);

                            numberWithPID++;
                            cellValuesVector.push_back(p);
                       }
                    }
                    if (!numberWithPID)
                    {
                        Info<<"foamToNASTRAN::main(): No faces found with PID: "
                            << pid <<" on surface: "
                            <<targetFileNames[surfI].lessExt()<<endl;
                    }
                    forceNumber++;
                   //targetFileNames[surfI]
                } // forAll(targetOutputPIDs)

        hermitianSpectra+= computeHermitianSpectra(allPointsGlobal, facesGlobal) ;
#ifdef DISPLAYERR
        std::vector<scalar> cellEtalonValues;
        Foam::readPrevPressureComputeError
        (
            "initialPressureData",
            cellEtalonValues
        );

        if (cellEtalonValues.size()!=cellValuesVector.size()) return 1;

        std::vector<scalar> cellProjectionErorrValues;
        scalar averageError(0.0);
        for (unsigned int ii=0; ii<cellEtalonValues.size();ii++)
        {
            scalar delta=( ( fabs(cellValuesVector[ii]-cellEtalonValues[ii]) )
                           /fabs( cellEtalonValues[ii] ) )*100.0;
            cellProjectionErorrValues.push_back(delta);
            averageError+=delta;
        }

        averageError/=cellEtalonValues.size();

        Info<<"foamToNASTRAN::main(): average projection error ="
            <<averageError<<" per cent."<<endl;

        for (unsigned int ii=0; ii<cellEtalonValues.size();ii++)
        {
            scalar err( cellProjectionErorrValues[ii] );
            if (err > 100.0 * averageError)
                cellProjectionErorrValues[ii] = averageError;
        }

        word file_name( targetFileNames[surfI] );
        if (isSourceSurfacesPresent) file_name=sourceFileNames[surfI];

        Foam::writeVTKModified
        (
            targetFileNames[surfI],
            allPointsGlobal,
            facesGlobal,
            cellProjectionErorrValues,
            false
         );
#endif

#ifdef DISPLAYDIST
                std::vector<scalar> nodalDistance;
        readPrevPointsComputeDist
        (
            "initialPointsData",
            allPointsGlobal,
            nodalDistance
        );

        word file_name( targetFileNames[surfI] );
        if (isSourceSurfacesPresent) file_name=sourceFileNames[surfI];

        Foam::writeVTKModified
        (
            file_name,
            allPointsGlobal,
            facesGlobal,
            nodalDistance,
            true
          );
#endif

            } //master && !parRun
            Info<<endl<<" foamToNASTRAN::main(): For surface="<<surfI
                <<" hermitian spectra="<<hermitianSpectra<<endl;
        } //forAll(targetFileNames, surfI)
        if (nFailedAttempts)
            Info<<endl<<"WARNING: foamToNASTRAN::main(): Detected "
                <<nFailedAttempts<<" failed projections."<<endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
