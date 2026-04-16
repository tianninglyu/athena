//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mhdturb_jeans.cpp
//! \brief Problem generator for configurable MHD turbulence with Jeans AMR criterion

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
Real density0;
Real pressure0;
Real vx0, vy0, vz0;
Real bx0, by0, bz0;

Real four_pi_g0 = 1.0;
Real njeans_refine = 4.0;
}  // namespace

//========================================================================================
//! \fn int JeansRefinementCondition(MeshBlock *pmb)
//! \brief Refine/derefine based on the minimum cell Jeans number in a MeshBlock
//========================================================================================
int JeansRefinementCondition(MeshBlock *pmb) {
  Real max_ratio = 0.0;

  // Assumes uniform cubic cells.
  const Real dx = pmb->pcoord->dx1f(0);
  const Real gconst = four_pi_g0/(4.0*PI);
  const Real jeans_j = 1.0/njeans_refine;
  const Real rho_fac = SQR(jeans_j)*PI/(gconst*SQR(dx));

  for (int k = pmb->ks-NGHOST; k <= pmb->ke+NGHOST; ++k) {
    for (int j = pmb->js-NGHOST; j <= pmb->je+NGHOST; ++j) {
      for (int i = pmb->is-NGHOST; i <= pmb->ie+NGHOST; ++i) {
        const Real rho = pmb->phydro->w(IDN, k, j, i);

        Real cs = 0.0;
        if (NON_BAROTROPIC_EOS) {
          const Real gamma = pmb->peos->GetGamma();
          cs = std::sqrt(gamma*pmb->phydro->w(IPR, k, j, i)/rho);
        } else {
          cs = pmb->peos->GetIsoSoundSpeed();
        }

        Real inv_beta = 0.0;
        if (MAGNETIC_FIELDS_ENABLED) {
          const Real b2 = SQR(pmb->pfield->bcc(IB1, k, j, i))
                        + SQR(pmb->pfield->bcc(IB2, k, j, i))
                        + SQR(pmb->pfield->bcc(IB3, k, j, i));
          if (b2 > 0.0) {
            Real pgas = 0.0;
            if (NON_BAROTROPIC_EOS) {
              pgas = pmb->phydro->w(IPR, k, j, i);
            } else {
              pgas = rho*SQR(cs);
            }
            const Real beta = 2.0*pgas/b2;
            if (beta > 0.0) inv_beta = 1.0/beta;
          }
        }

        const Real ceff = cs*std::sqrt(1.0 + 0.74*inv_beta);
        const Real rho_crit = rho_fac*ceff;
        const Real ratio = rho/rho_crit;
        max_ratio = std::max(max_ratio, ratio);
      }
    }
  }

  if (max_ratio > 1.0) return 1;
  if (max_ratio < 1.0/2.5) return -1;
  return 0;
}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Initialize problem-specific data that need to be shared across mesh blocks
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (SELF_GRAVITY_ENABLED) {
    four_pi_g0 = pin->GetReal("problem", "four_pi_G");
    SetFourPiG(four_pi_g0);

    if (adaptive) {
      njeans_refine = pin->GetOrAddReal("problem", "njeans", 4.0);
      EnrollUserRefinementCondition(JeansRefinementCondition);
    }
  }

  density0 = pin->GetReal("problem", "density");
  pressure0 = pin->GetReal("problem", "pressure");
  vx0 = pin->GetOrAddReal("problem", "vx", 0.0);
  vy0 = pin->GetOrAddReal("problem", "vy", 0.0);
  vz0 = pin->GetOrAddReal("problem", "vz", 0.0);

  if (MAGNETIC_FIELDS_ENABLED) {
    bx0 = pin->GetOrAddReal("problem", "bx", 0.0);
    by0 = pin->GetOrAddReal("problem", "by", 0.0);
    bz0 = pin->GetOrAddReal("problem", "bz", 0.0);
  } else {
    bx0 = 0.0;
    by0 = 0.0;
    bz0 = 0.0;
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initialize conserved variables for the MHD turbulence problem
//========================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gm1 = 0.0;
  if (NON_BAROTROPIC_EOS) {
    gm1 = peos->GetGamma() - 1.0;
  }

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = density0;
        phydro->u(IM1,k,j,i) = density0*vx0;
        phydro->u(IM2,k,j,i) = density0*vy0;
        phydro->u(IM3,k,j,i) = density0*vz0;

        if (NON_BAROTROPIC_EOS) {
          Real kinetic = 0.5*density0*(SQR(vx0) + SQR(vy0) + SQR(vz0));
          phydro->u(IEN,k,j,i) = pressure0/gm1 + kinetic;
        }
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = bx0;
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = by0;
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = bz0;
        }
      }
    }

    if (NON_BAROTROPIC_EOS) {
      Real magnetic = 0.5*(SQR(bx0) + SQR(by0) + SQR(bz0));
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            phydro->u(IEN,k,j,i) += magnetic;
          }
        }
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief No post-processing required for this problem
//========================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}
