//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mhdturb.cpp
//! \brief Problem generator for configurable MHD turbulence setups

// C headers

// C++ headers
#include <cmath>

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
}  // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Initialize problem-specific data that need to be shared across mesh blocks
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (SELF_GRAVITY_ENABLED) {
    Real four_pi_G = pin->GetReal("problem", "four_pi_G");
    SetFourPiG(four_pi_G);
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

