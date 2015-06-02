// General relativistic Fishbone-Moncrief torus generator

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <cmath>      // exp(), pow(), sin(), sqrt()
#include <algorithm>  // max()

// Athena headers
#include "../athena.hpp"                   // enums, Real
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../bvals/bvals.hpp"              // BoundaryValues
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../fluid/fluid.hpp"              // Fluid
#include "../fluid/eos/eos.hpp"            // FluidEqnOfState

// Declarations
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
void FixedTop(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
void FixedBottom(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke);
static void reset_l_from_r_peak();
static Real log_h_aux(Real r, Real sin_theta);
static void calculate_velocity(Real r, Real sin_theta, Real *pu0, Real *pu3);
static void set_state(
    Real rho, Real pgas, Real v1, Real v2, Real v3, int k, int j, int i,
    AthenaArray<Real> &prim, AthenaArray<Real> &prim_half);
static void set_conserved_cell(MeshBlock *pmb, Real gamma_adi, int k, int j, int i,
    AthenaArray<Real> &cons);

// Global variables
static Real m, a;                            // black hole parameters
static Real gamma_adi, k_adi;                // fluid parameters
static Real r_edge, r_peak, l, rho_max;      // disk parameters
static Real rho_min, rho_pow, u_min, u_pow;  // background parameters

// Function for setting initial conditions
// Inputs:
//   pfl: Fluid
//   pfd: Field (unused)
//   pin: parameters
// Outputs: (none)
// Notes:
//   initializes Fishbone-Moncrief torus
//     sets both primitive and conserved variables
//   defines and enrolls fixed r- and theta-direction boundary conditions
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
//              Fishbone 1977, ApJ 215 323 (F)
//   assumes x3 is axisymmetric direction
// TODO: only works in Schwarzschild (assumed metric)
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  // Prepare index bounds
  MeshBlock *pb = pfl->pmy_block;
  int il = pb->is - NGHOST;
  int iu = pb->ie + NGHOST;
  int jl = pb->js;
  int ju = pb->je;
  if (pb->block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = pb->ks;
  int ku = pb->ke;
  if (pb->block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Get mass and spin of black hole
  m = pb->pcoord->GetMass();
  a = pb->pcoord->GetSpin();

  // Get ratio of specific heats
  gamma_adi = pfl->pf_eos->GetGamma();

  // Read other properties
  k_adi = pin->GetReal("problem", "k_adi");
  r_edge = pin->GetReal("problem", "r_edge");
  r_peak = pin->GetReal("problem", "r_peak");
  l = pin->GetReal("problem", "l");
  rho_max = pin->GetReal("problem", "rho_max");
  rho_min = pin->GetReal("problem", "rho_min");
  rho_pow = pin->GetReal("problem", "rho_pow");
  u_min = pin->GetReal("problem", "u_min");
  u_pow = pin->GetReal("problem", "u_pow");
  potential_cutoff = pin->GetReal("problem", "cutoff");

  // Reset l if valid r_peak given
  if (r_peak >= 0.0)
    reset_l_from_r_peak();

  // Prepare to keep track of which cells are in the torus for normalization purposes
  Real rho_peak = 0.0;
  AthenaArray<bool> in_torus;
  in_torus.NewAthenaArray(ju+1, iu+1);

  // Initialize primitive values
  Real log_h_edge = log_h_aux(r_edge, 1.0);
  for (int j = jl; j <= ju; j++)
    for (int i = il; i <= iu; i++)
    {
      // Get Boyer-Lindquist coordinates of cell
      Real r, theta, phi;
      pmb->pcoord->GetBoyerLindquistCoordinates(pb->x1v(i), pb->x2v(j), pb->x3v(kl),
          &r, &theta, &phi);
      Real sin_theta = std::sin(theta);

      // Determine if we are in the torus
      Real log_h;
      in_torus(j,i) = false;
      if (r >= r_edge)
      {
        log_h = log_h_aux(r, sin_theta) - log_h_edge;  // (FM 3.6)
        if (log_h >= 0.0)
          in_torus(j,i) = true;
      }

      // Calculate thermodynamic quantities
      Real rho, pgas;
      if (in_torus(j,i))
      {
        Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1);
        rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0));
        pgas = pgas_over_rho * rho;
        rho_peak = std::max(rho_peak, rho);
      }
      else
      {
        //rho = rho_min * std::pow(r/r_edge, rho_pow);
        //Real u = u_min * std::pow(r/r_edge, u_pow);
        rho = rho_min;
        Real u = u_min;
        pgas = (gamma_adi-1.0) * u;
      }

      // Calculate velocity
      Real u0, u3;
      calculate_velocity(r, sin_theta, &u0, &u3);

      // Convert velocity back to preferred coordinate system
      Real u0_pref, u1_pref, u2_pref, u3_pref;
      pb->pcoord->TransformVectorCell(u0, 0.0, 0.0, u3, kl, j, i,
          &u0_pref, &u1_pref, &u2_pref, &u3_pref);
      Real v3 = u3_pref/u0_pref;

      // Set primitive values
      for (int k = kl; k <= ku; k++)
        set_state(rho, pgas, 0.0, 0.0, v3, k, j, i, pfl->w, pfl->w1);
    }

  // Initialize magnetic fields
  if (MAGNETIC_FIELDS_ENABLED)
  {
    // Prepare 2D arrays of vector potential values
    AthenaArray<Real> a_phi_cells, a_phi_edges;
    a_phi_cells.NewAthenaArray(ju+1, iu+1);
    a_phi_edges.NewAthenaArray(ju+2, iu+2);

    // Go through 2D slice, setting vector potential in cells
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        // Get Boyer-Lindquist coordinates
        Real r, theta, phi;
        pmb->pcoord->GetBoyerLindquistCoordinates(pb->x1v(i), pb->x2v(j), pb->x3v(kl),
            &r, &theta, &phi);

        // Calculate A_phi as proportional to rho
        if (r >= r_edge)
        {
          Real log_h = log_h_aux(r, sin_theta) - log_h_edge;  // (FM 3.6)
          if (log_h >= 0.0)
          {
            Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1);
            Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0));
            rho_peak = std::max(rho_peak, rho);
            a_phi_cells(j,i) = rho;
          }
        }
      }

    // Go through 2D slice, setting vector potential at edges
    for (int j = jl; j <= ju+1; j++)
      for (int i = il; i <= iu+1; i++)
      {
        // Get Boyer-Lindquist coordinates
        Real r, theta, phi;
        pmb->pcoord->GetBoyerLindquistCoordinates(pb->x1f(i), pb->x2f(j), pb->x3v(kl),
            &r, &theta, &phi);

        // Calculate A_phi as proportional to rho
        if (r >= r_edge)
        {
          Real log_h = log_h_aux(r, sin_theta) - log_h_edge;  // (FM 3.6)
          if (log_h >= 0.0)
          {
            Real pgas_over_rho = (gamma_adi-1.0)/gamma_adi * (std::exp(log_h)-1);
            Real rho = std::pow(pgas_over_rho/k_adi, 1.0/(gamma_adi-1.0));
            rho_peak = std::max(rho_peak, rho);
            a_phi_edges(j,i) = rho;
          }
        }
      }

    // Normalize vector potential
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
        a_phi_cells(j,i) = std::max(a_phi_cells(j,i)/rho_peak - potential_cutoff, 0.0);
    for (int j = jl; j <= ju+1; j++)
      for (int i = il; i <= iu+1; i++)
        a_phi_edges(j,i) = std::max(a_phi_edges(j,i)/rho_peak - potential_cutoff, 0.0);

    // Set magnetic fields according to vector potential
    // Note: This does very rough differencing for in-face fields on exterior faces of
    //    domain. This should not matter, as these will be identically 0 in nice
    //    coordinate systems or as long as the initial torus is within the domain.
    for (int k = kl; k <= ku+1; ++k)
      for (int j = jl; j <= ju+1; ++j)
        for (int i = il; i <= iu+1; ++i)
        {
          // Set B^1
          if (j != ju+1 and k != ku+1)
          {
            Real r, theta, phi;
            pmb->pcoord->GetBoyerLindquistCoordinates(
                pb->x1f(i), pb->x2v(j), pb->x3v(k), &r, &theta, &phi);
            Real r_1, theta_1, phi_1;
            pmb->pcoord->GetBoyerLindquistCoordinates(
                pb->x1f(i), pb->x2f(j), pb->x3v(k), &r_1, &theta_1, &phi_1);
            Real r_2, theta_2, phi_2;
            pmb->pcoord->GetBoyerLindquistCoordinates(
                pb->x1f(i), pb->x2f(j+1), pb->x3v(k), &r_2, &theta_2, &phi_2);
            Real bbr = (a_phi_edges(j+1,i) - a_phi_edges(j,i)) / (theta_2 - theta_1);
            Real a_phi_1, a_phi_2;
            if (i == il)
            {
              a_phi_1 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j+1,i));
              a_phi_2 = a_phi_cells(j,i);
              r_1 = r;
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i), pb->x2v(j), pb->x3v(k), &r_2, &theta_2, &phi_2);
            }
            else if (i == iu+1)
            {
              a_phi_1 = a_phi_cells(j,i-1);
              a_phi_2 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j+1,i));
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i-1), pb->x2v(j), pb->x3v(k), &r_1, &theta_1, &phi_1);
              r_2 = r;
            }
            else
            {
              a_phi_1 = a_phi_cells(j,i-1);
              a_phi_2 = a_phi_cells(j,i);
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i-1), pb->x2v(j), pb->x3v(k), &r_1, &theta_1, &phi_1);
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i), pb->x2v(j), pb->x3v(k), &r_2, &theta_2, &phi_2);
            }
            Real bbtheta = (a_phi_2 - a_phi_1) / (r_2 - r_1);
            Real ut, uphi;
            Real sin_theta = std::sin(theta);
            calculate_velocity(r, sin_theta, &ut, &uphi);
            Real sin_sq_theta = SQR(sin_theta);
            Real cos_sq_theta = 1.0 - sin_sq_theta;
            Real rho_sq = SQR(r) + SQR(a) * cos_sq_theta;
            bt = -2.0*m*a*r * SQR(sin_theta) / rho_sq * bbr * ut;
            br = 1.0/ut * bbr;
            btheta = 1.0/ut * bbtheta;
            Real u0, u1, u2, u3;
            pmb->pcoord->TransformVectorFace1(ut, 0.0, 0.0, uphi, k, j, i,
                &u0, &u1, &u2, &u3);
            Real b0, b1, b2, b3;
            pmb->pcoord->TransformVectorFace1(bt, br, btheta, 0.0, k, j, i,
                &b0, &b1, &b2, &b3);
            pfd->b.x1f(k,j,i) = b1 * u0 - b0 * u1;
          }

          // Set B^2
          if (i != iu+1 and k != ku+1)
          {
            Real r, theta, phi;
            pmb->pcoord->GetBoyerLindquistCoordinates(
                pb->x1v(i), pb->x2f(j), pb->x3v(k), &r, &theta, &phi);
            Real r_1, theta_1, phi_1;
            pmb->pcoord->GetBoyerLindquistCoordinates(
                pb->x1f(i), pb->x2f(j), pb->x3v(k), &r_1, &theta_1, &phi_1);
            Real r_2, theta_2, phi_2;
            pmb->pcoord->GetBoyerLindquistCoordinates(
                pb->x1f(i+1), pb->x2f(j), pb->x3v(k), &r_2, &theta_2, &phi_2);
            Real bbtheta = (a_phi_edges(j,i+1) - a_phi_edges(j,i)) / (r_2 - r_1);
            Real a_phi_1, a_phi_2;
            if (j == jl)
            {
              a_phi_1 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j,i+1));
              a_phi_2 = a_phi_cells(j,i);
              theta_1 = theta;
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i), pb->x2v(j), pb->x3v(k), &r_2, &theta_2, &phi_2);
            }
            else if (j == ju+1)
            {
              a_phi_1 = a_phi_cells(j-1,i);
              a_phi_2 = 0.5 * (a_phi_edges(j,i) + a_phi_edges(j,i+1));
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i), pb->x2v(j-1), pb->x3v(k), &r_1, &theta_1, &phi_1);
              theta_2 = theta;
            }
            else
            {
              a_phi_1 = a_phi_cells(j-1,i);
              a_phi_2 = a_phi_cells(j,i);
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i), pb->x2v(j-1), pb->x3v(k), &r_1, &theta_1, &phi_1);
              pmb->pcoord->GetBoyerLindquistCoordinates(
                  pb->x1v(i), pb->x2v(j), pb->x3v(k), &r_2, &theta_2, &phi_2);
            }
            Real bbr = (a_phi_2 - a_phi_1) / (theta_2 - theta_1);
            Real ut, uphi;
            Real sin_theta = std::sin(theta);
            calculate_velocity(r, sin_theta, &ut, &uphi);
            Real sin_sq_theta = SQR(sin_theta);
            Real cos_sq_theta = 1.0 - sin_sq_theta;
            Real rho_sq = SQR(r) + SQR(a) * cos_sq_theta;
            bt = -2.0*m*a*r * SQR(sin_theta) / rho_sq * bbr * ut;
            br = 1.0/ut * bbr;
            btheta = 1.0/ut * bbtheta;
            Real u0, u1, u2, u3;
            pmb->pcoord->TransformVectorFace2(ut, 0.0, 0.0, uphi, k, j, i,
                &u0, &u1, &u2, &u3);
            Real b0, b1, b2, b3;
            pmb->pcoord->TransformVectorFace2(bt, br, btheta, 0.0, k, j, i,
                &b0, &b1, &b2, &b3);
            pfd->b.x2f(k,j,i) = b2 * u0 - b0 * u2;
          }

          // Set B^3
          if (i != iu+1 and j != ju+1)
            pfd->b.x3f(k,j,i) = 0.0;
        }
    a_phi_cells.DeleteAthenaArray();
    a_phi_edges.DeleteAthenaArray();
  }

  // Normalize density and pressure
  if (rho_max > 0.0)
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i)
          if (in_torus(j,i))
          {
            pfl->w(IDN,k,j,i) /= rho_peak;
            pfl->w1(IDN,k,j,i) /= rho_peak;
            pfl->w(IEN,k,j,i) /= rho_peak;
            pfl->w1(IEN,k,j,i) /= rho_peak;
            // TODO: better way to induce splash and crash
            //w(IM3,k,j,i) /= 2.0;
          }
  in_torus.DeleteAthenaArray();

  // Initialize conserved values
  // TODO: allow for magnetic fields to be nonzero
  AthenaArray<Real> b;
  b.NewAthenaArray(ku, ju, iu);
  pb->pcoord->PrimToCons(pfl->w, b, gamma_adi/(gamma_adi-1.0), pfl->u);
  b.DeleteAthenaArray();

  // Enroll boundary functions
  // TODO enroll field boundary conditions
  pb->pbval->EnrollFluidBoundaryFunction(inner_x1, FixedInner);
  pb->pbval->EnrollFluidBoundaryFunction(outer_x1, FixedOuter);
  pb->pbval->EnrollFluidBoundaryFunction(inner_x2, FixedTop);
  pb->pbval->EnrollFluidBoundaryFunction(outer_x2, FixedBottom);
  return;
}

// Inner boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along inner x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedInner(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Extract boundary indices
  int il = is - NGHOST;
  int iu = is;
  int jl = js;
  int ju = je;
  int kl = ks;
  int ku = ke;

  // Set conserved values
  Real r = pmb->x1v(iu);
  Real d, e;
  calculate_conserved(r, d, e);
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Outer boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along outer x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedOuter(MeshBlock *pmb, AthenaArray<Real> &cons,
                int is, int ie, int js, int je, int ks, int ke)
{
  // Extract boundary indices
  int il = ie;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  int kl = ks;
  int ku = ke;

  // Set conserved values
  Real r = pmb->x1v(il);
  Real d, e;
  calculate_conserved(r, d, e);
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Top boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along inner x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedTop(MeshBlock *pmb, AthenaArray<Real> &cons,
              int is, int ie, int js, int je, int ks, int ke)
{
  // Extract boundary indices
  int il = is;
  int iu = ie;
  int jl = js - NGHOST;
  int ju = js;
  int kl = ks;
  int ku = ke;

  // Set conserved values
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        Real r = pmb->x1v(i);
        Real d, e;
        calculate_conserved(r, d, e);
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Bottom boundary condition
// Inputs:
//   pmb: pointer to block
// Outputs:
//   cons: conserved quantities set along inner x1-boundary
// Notes:
//   references Hawley, Smarr, & Wilson 1984, ApJ 277 296 (HSW)
// TODO: only works in Schwarzschild (assumed metric)
void FixedBottom(MeshBlock *pmb, AthenaArray<Real> &cons,
                 int is, int ie, int js, int je, int ks, int ke)
{
  // Extract boundary indices
  int il = is;
  int iu = ie;
  int jl = je;
  int ju = je + NGHOST;
  int kl = ks;
  int ku = ke;

  // Set conserved values
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
      for (int i = il; i <= iu; i++)
      {
        Real r = pmb->x1v(i);
        Real d, e;
        calculate_conserved(r, d, e);
        cons(IDN,k,j,i) = d;
        cons(IEN,k,j,i) = e;
        cons(IM1,k,j,i) = 0.0;
        cons(IM2,k,j,i) = 0.0;
        cons(IM3,k,j,i) = 0.0;
      }
  return;
}

// Function for calculating angular momentum variable l
// Inputs: (none)
// Outputs:
//   sets global variable l = u^t u_\phi such that pressure maximum occurs at r_peak
// Notes:
//   beware many different definitions of l abound
//     this is *not* -u_phi/u_t
//   Harm has a similar function: lfish_calc() in init.c
//     Harm's function assumes M = 1 and that corotation is desired
//     it is equivalent to this, though seeing this requires much manipulation
//   implements (3.8) from Fishbone & Moncrief 1976, ApJ 207 962
//   assumes corotation
//   TODO: add counterrotation option
static void reset_l_from_r_peak()
{
  Real num = SQR(SQR(r_peak)) + SQR(a*r_peak) - 2.0*m*SQR(a)*r_peak
      - a*(SQR(r_peak)-SQR(a))*std::sqrt(m*r_peak);
  Real denom = SQR(r_peak) - 3.0*m*r_peak + 2.0*a*std::sqrt(m*r_peak);
  l = 1.0/r_peak * std::sqrt(m/r_peak) * num/denom;
  return;
}

// Function for helping to calculate enthalpy
// Inputs:
//   r: radial Boyer-Lindquist coordinate
//   sin_sq_theta: square of sine of polar Boyer-Lindquist coordinate
// Outputs:
//   returned value: log(h)
// Notes:
//   enthalpy defined here as h = p_gas/rho
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
//   implements first half of (FM 3.6)
static Real log_h_aux(Real r, Real sin_theta)
{
  Real sin_sq_theta = SQR(sin_theta);
  Real cos_sq_theta = 1.0 - sin_sq_theta;
  Real delta = SQR(r) - 2.0*m*r + SQR(a);                // \Delta
  Real sigma = SQR(r) + SQR(a)*cos_sq_theta;             // \Sigma
  Real aa = SQR(SQR(r)+SQR(a))
      - delta*SQR(a)*sin_sq_theta;                       // A
  Real exp_2nu = sigma * delta / aa;                     // \exp(2\nu) (FM 3.5)
  Real exp_2psi = aa / sigma * sin_sq_theta;             // \exp(2\psi) (FM 3.5)
  Real exp_neg2chi = exp_2nu / exp_2psi;                 // \exp(-2\chi) (cf. FM 2.15)
  Real omega = 2.0*m*a*r/aa;                             // \omega (FM 3.5)
  Real var_a = std::sqrt(1.0 + 4.0*SQR(l)*exp_neg2chi);
  Real var_b = 0.5 * std::log((1.0+var_a)
      / (sigma*delta/aa));
  Real var_c = -0.5 * var_a;
  Real var_d = -l * omega;
  return var_b + var_c + var_d;                          // (FM 3.4)
}

// Function for computing 4-velocity components at a given position
// Inputs:
//   r: Boyer-Lindquist r
//   sin_theta: sine of Boyer-Lindquist theta
// Outputs:
//   pu0: u^t set (Boyer-Lindquist coordinates)
//   pu3: u^\phi set (Boyer-Lindquist coordinates)
// Notes:
//   The formula for u^3 as a function of u_{(\phi)} is tedious to derive,
//       but this matches the formula used in Harm (init.c).
static void calculate_velocity(Real r, Real sin_theta, Real *pu0, Real *pu3)
{
  Real sin_sq_theta = SQR(sin_theta);
  Real cos_sq_theta = 1.0 - sin_sq_theta;
  Real delta = SQR(r) - 2.0*m*r + SQR(a);            // \Delta
  Real sigma = SQR(r) + SQR(a)*cos_sq_theta;         // \Sigma
  Real aa = SQR(SQR(r)+SQR(a))
      - delta*SQR(a)*sin_sq_theta;                   // A
  Real exp_2nu = sigma * delta / aa;                 // \exp(2\nu) (FM 3.5)
  Real exp_2psi = aa / sigma * sin_sq_theta;         // \exp(2\psi) (FM 3.5)
  Real exp_neg2chi = exp_2nu / exp_2psi;             // \exp(-2\chi) (cf. FM 2.15)
  Real u_phi_proj_a = 1.0 + 4.0*SQR(l)*exp_neg2chi;
  Real u_phi_proj_b = -1.0
      + std::sqrt(u_phi_proj_a);
  Real u_phi_proj = std::sqrt(0.5 * u_phi_proj_b);   // (FM 3.3)
  Real u_3 = std::sqrt(aa/sigma) * sin_theta
      * u_phi_proj;                                  // (FM 2.12, F 2.5, FM 3.5)
  Real u3_a = (1.0+SQR(u_phi_proj))
      / (aa*sigma*delta);
  Real u3_b = 2.0*m*a*r * std::sqrt(u3_a);
  Real u3_c = std::sqrt(sigma/aa) / sin_theta;
  Real u3 = u3_b + u3_c * u_phi_proj;
  Real g_00 = -(1.0 - 2.0*m*r/sigma);
  Real g_03 = -2.0*m*a*r/sigma * sin_sq_theta;
  Real g_33 = (sigma + (1.0 + 2.0*m*r/sigma)
      * SQR(a) * sin_sq_theta) * sin_sq_theta;
  Real u0_a = (SQR(g_03) - g_00*g_33) * SQR(u3);
  Real u0_b = std::sqrt(u0_a - g_00);
  Real u0 = -1.0/g_00 * (g_03*u3 + u0_b);
  *pu0 = u0;
  *pu3 = u3;
  return;
}

// Function for setting conserved variables in a cell given the primitives
// Inputs:
//   rho: density
//   pgas: gas pressure
//   v1,v2,v3: 3-velocity components
//   k,j,i: indices for cell to be set
// Outputs:
//   prim, prim_half: primitive values set
static void set_state(
    Real rho, Real pgas, Real v1, Real v2, Real v3, int k, int j, int i,
    AthenaArray<Real> &prim, AthenaArray<Real> &prim_half)
{
  prim(IDN,k,j,i) = prim_half(IDN,k,j,i) = rho;
  prim(IEN,k,j,i) = prim_half(IEN,k,j,i) = pgas;
  prim(IM1,k,j,i) = prim_half(IM1,k,j,i) = v1;
  prim(IM2,k,j,i) = prim_half(IM2,k,j,i) = v2;
  prim(IM3,k,j,i) = prim_half(IM3,k,j,i) = v3;
  return;
}

// Function for setting boundary conserved quantities in a given cell
// Inputs:
//   pmb: pointer to block
//   gamma_adi: ratio of specific heats \Gamma
//   k,j,i: indices of cell in which conserved quantities are to be set
// Outputs:
//   cons: conserved quantities set in desired cell
// Notes:
//   first implements same procedure as ProblemGenerator()
//   then converts primitive to conserved values
//   references Fishbone & Moncrief 1976, ApJ 207 962 (FM)
//              Fishbone 1977, ApJ 215 323 (F)
//   TODO: update for MHD
//   TODO: investigate way of storing desired values rather than recomputing them
static void set_conserved_cell(MeshBlock *pmb, Real gamma_adi, int k, int j, int i,
    AthenaArray<Real> &cons)
{
  // Get Boyer-Lindquist coordinates of cell
  Real r, theta, phi;
  pmb->pcoord->GetBoyerLindquistCoordinates(pmb->x1v(i), pmb->x2v(j), pmb->x3v(k),
      &r, &theta, &phi);
  Real sin_theta = std::sin(theta);
  Real sin_sq_theta = SQR(sin_theta);
  Real cos_sq_theta = 1.0 - sin_sq_theta;

  // Calculate thermodynamic quantities
  //Real rho = rho_min * std::pow(r/r_edge, rho_pow);
  //Real u = u_min * std::pow(r/r_edge, u_pow);
  Real rho = rho_min;
  Real u = u_min;
  Real pgas = (gamma_adi-1.0) * u;

  // Calculate velocity
  Real delta = SQR(r) - 2.0*m*r + SQR(a);            // \Delta
  Real sigma = SQR(r) + SQR(a)*cos_sq_theta;         // \Sigma
  Real aa = SQR(SQR(r)+SQR(a))
      - delta*SQR(a)*sin_sq_theta;                   // A
  Real exp_2nu = sigma * delta / aa;                 // \exp(2\nu) (FM 3.5)
  Real exp_2psi = aa / sigma * sin_sq_theta;         // \exp(2\psi) (FM 3.5)
  Real exp_neg2chi = exp_2nu / exp_2psi;             // \exp(-2\chi) (cf. FM 2.15)
  Real u_phi_proj_a = 1.0 + 4.0*SQR(l)*exp_neg2chi;
  Real u_phi_proj_b = -1.0
      + std::sqrt(u_phi_proj_a);
  Real u_phi_proj = std::sqrt(0.5 * u_phi_proj_b);   // (FM 3.3)
  Real u_3 = std::sqrt(aa/sigma) * sin_theta
      * u_phi_proj;                                  // (FM 2.12, F 2.5, FM 3.5)
  Real u3_a = (1.0+SQR(u_phi_proj))
      / (aa*sigma*delta);
  Real u3_b = 2.0*m*a*r * std::sqrt(u3_a);
  Real u3_c = std::sqrt(sigma/aa) / sin_theta;
  Real u3 = u3_b + u3_c * u_phi_proj;
  Real g_00 = -(1.0 - 2.0*m*r/sigma);
  Real g_03 = -2.0*m*a*r/sigma * sin_sq_theta;
  Real g_33 = (sigma + (1.0 + 2.0*m*r/sigma)
      * SQR(a) * sin_sq_theta) * sin_sq_theta;
  Real u0_a = (SQR(g_03) - g_00*g_33) * SQR(u3);
  Real u0_b = std::sqrt(u0_a - g_00);
  Real u0 = -1.0/g_00 * (g_03*u3 + u0_b);

  // Convert velocity back to preferred coordinate system
  Real u0_pref, u1_pref, u2_pref, u3_pref;
  Real u_0_pref, u_1_pref, u_2_pref, u_3_pref;
  pmb->pcoord->TransformVectorCell(u0, 0.0, 0.0, u3, k, j, i,
      &u0_pref, &u1_pref, &u2_pref, &u3_pref);
  pmb->pcoord->LowerVectorCell(u0_pref, u1_pref, u2_pref, u3_pref, k, j, i,
      &u_0_pref, &u_1_pref, &u_2_pref, &u_3_pref);

  // Set conserved values
  Real gamma_adi_red = gamma_adi/(gamma_adi-1.0);
  Real w = rho + gamma_adi_red * pgas;
  Real ptot = pgas;
  cons(IDN,k,j,i) = rho * u0_pref;
  cons(IEN,k,j,i) = w * u0_pref * u_0_pref + ptot;
  cons(IM1,k,j,i) = w * u0_pref * u_1_pref;
  cons(IM2,k,j,i) = w * u0_pref * u_2_pref;
  cons(IM3,k,j,i) = w * u0_pref * u_3_pref;
  return;
}
