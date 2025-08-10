#include "checkMeshQuality.H"
#include "meshes/polyMesh/polyMesh.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "motionSmoother/motionSmoother.H"
#include "sampledSurface/writers/surfaceWriter.H"
#include "checkTools.H"
#include "motionSmoother/polyMeshGeometry/meshStatistics.H"

Foam::label Foam::checkMeshQuality
(
    const fvMesh& mesh,
    const dictionary& dict,
    const wordList& metrics,
    const bool& report,
    const bool& writeCombinedMetric,
    const autoPtr<surfaceWriter>& writer
)
{
    label noFailedChecks = 0;

    {
        PtrList<meshStatistics> metricsFields(0);

        if (metrics.size())
        {
            metricsFields.setSize(metrics.size());
            forAll(metrics, metricI)
            {
                word metricName = metrics[metricI];
                dictionary histSettings;
                if (dict.found("histograms"))
                {
                    dictionary allHist = dict.subDict("histograms");
                    if (allHist.found(metricName))
                    {
                        histSettings = allHist.subDict(metricName);
                    }
                }

                metricsFields.set
                (
                    metricI,
                    new meshStatistics
                    (
                        mesh,
                        metricName,
                        histSettings
                    )
                );
            }
            if (writeCombinedMetric)
            {
                metricsFields.setSize(metrics.size()+1);
                word metricName("meshQuality");
                dictionary histSettings;
                if (dict.found("histograms"))
                {
                    dictionary allHist = dict.subDict("histograms");
                    if (allHist.found(metricName))
                    {
                        histSettings = allHist.subDict(metricName);
                    }
                }

                metricsFields.set
                (
                    metrics.size(),
                    new meshStatistics
                    (
                        mesh,
                        metricName,
                        histSettings
                    )
                );
            }
        }

        faceSet faces(mesh, "meshQualityFaces", mesh.nFaces()/100+1);
        motionSmoother::checkMesh
        (
            false,
            mesh,
            dict,
            faces,
            report,
            &metricsFields
        );

        label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            noFailedChecks++;

            Info<< "  <<Writing " << nFaces
                << " faces in error to set " << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (writer.valid())
            {
                mergeAndWrite(writer(), faces);
            }
        }

        if (metrics.size())
        {
            Info<<"\nOutputting mesh quality fields... "<<nl<<endl;

            forAll(metrics, metricI)
            {
                Info<<"Writing field: "<<metrics[metricI]<<endl;
                metricsFields[metricI].field().write();
                metricsFields[metricI].writeStats();
            }

            if (writeCombinedMetric)
            {
                volScalarField& mf = metricsFields[metrics.size()].field();
                Info<<"Writing field: "<<mf.name()<<endl;
                mf /= metrics.size();
                mf.write();
                metricsFields[metrics.size()].calcStats(mf);
                metricsFields[metrics.size()].writeStats();
            }
        }
    }

    return noFailedChecks;
}
