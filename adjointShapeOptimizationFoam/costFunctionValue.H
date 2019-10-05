
J = 0.0;
// integrate in the boundary patches
forAll(mesh.boundary(), patchi)
{
    // only over non-empty patches
    if (!p.boundaryField()[patchi].empty())
    {
        J -= gSum(
            (p.boundaryField()[patchi] + 0.5*pow(mag(U.boundaryField()[patchi]),2) )
            *( U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi])
        );
    }
}
