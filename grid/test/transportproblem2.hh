//===============================================================
// We want to solve the PDE
//
//   dc/dt + div( uc ) = 0   in Omega
//                   c = b   on Gamma_in
//                   c = c0  for t<0
//
// We define functions for the parameters c0, u and b.
//===============================================================

// the initial condition c0
template<int dimworld, class ct>
double c0 (const Dune::FieldVector<ct,dimworld>& x)
{
  if (x.two_norm()>0.125 && x.two_norm()<0.5)
	return 1.0;
  else
	return 0.0;
}

// the boundary condition b on inflow boundary
template<int dimworld, class ct>
double b (const Dune::FieldVector<ct,dimworld>& x, double t)
{
  return 0.0;
  if (x.two_norm()<t+0.125)
	return 1.0;
  else
	return 0.0;
}

// the vector field u is returned in r
template<int dimworld, class ct>
Dune::FieldVector<double,dimworld> u (const Dune::FieldVector<ct,dimworld>& x, double t)
{
  Dune::FieldVector<double,dimworld> r(1.0);
  return r;
}
