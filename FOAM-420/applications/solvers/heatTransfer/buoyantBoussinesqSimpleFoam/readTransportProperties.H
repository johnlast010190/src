singlePhaseTransportModel laminarTransport(U, phi);

// Reference temperature [K]
dimensionedScalar TRef
(
    "TRef",
    dimTemperature,
    gSum(T.primitiveField()*mesh.V())/gSum(mesh.V())
);

if (laminarTransport.found("TRef"))
{
    TRef = dimensionedScalar
    (
        "TRef", dimTemperature, laminarTransport.lookup("TRef")
    );
}

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    laminarTransport.found("beta") ?
    dimensionedScalar
    (
        "beta",
        dimless/dimTemperature,
        laminarTransport.lookup("beta")
    ) :
    1/TRef
);
