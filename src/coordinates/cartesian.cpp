//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Primary header
#include "coordinates.hpp"

// Athena++ headers
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh.hpp"            // MeshBlock

//======================================================================================
//! \file cartesian.cpp
//  \brief implements Coordinates class functions for Cartesian (x-y-z) coordinates
//======================================================================================

//--------------------------------------------------------------------------------------
// Coordinates constructor

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // initialize volume-averaged positions and spacing
  // x1-direction: x1v = dx/2
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pmb->x1v(i) = 0.5*(pmb->x1f(i+1) + pmb->x1f(i));
  }
  for (int i=is-(NGHOST); i<=ie+(NGHOST)-1; ++i) {
    pmb->dx1v(i) = pmb->x1v(i+1) - pmb->x1v(i);
  }

  // x2-direction: x2v = dy/2
  if (pmb->block_size.nx2 == 1) {
    pmb->x2v(js) = 0.5*(pmb->x2f(js+1) + pmb->x2f(js));
    pmb->dx2v(js) = pmb->dx2f(js);
  } else {
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
      pmb->x2v(j) = 0.5*(pmb->x2f(j+1) + pmb->x2f(j));
    }
    for (int j=js-(NGHOST); j<=je+(NGHOST)-1; ++j) {
      pmb->dx2v(j) = pmb->x2v(j+1) - pmb->x2v(j);
    }
  }

  // x3-direction: x3v = dz/2
  if (pmb->block_size.nx3 == 1) {
    pmb->x3v(ks) = 0.5*(pmb->x3f(ks+1) + pmb->x3f(ks));
    pmb->dx3v(ks) = pmb->dx3f(ks);
  } else {
    for (int k=ks-(NGHOST); k<=ke+(NGHOST); ++k) {
      pmb->x3v(k) = 0.5*(pmb->x3f(k+1) + pmb->x3f(k));
    }
    for (int k=ks-(NGHOST); k<=ke+(NGHOST)-1; ++k) {
      pmb->dx3v(k) = pmb->x3v(k+1) - pmb->x3v(k);
    }
  }

  if(pmb->pmy_mesh->multilevel==true) { // calc coarse coodinates
    int cis = pmb->cis; int cjs = pmb->cjs; int cks = pmb->cks;
    int cie = pmb->cie; int cje = pmb->cje; int cke = pmb->cke;
    for (int i=cis-(pmb->cnghost); i<=cie+(pmb->cnghost); ++i) {
      pmb->coarse_x1v(i) = 0.5*(pmb->coarse_x1f(i+1) + pmb->coarse_x1f(i));
    }
    for (int i=cis-(pmb->cnghost); i<=cie+(pmb->cnghost)-1; ++i) {
      pmb->coarse_dx1v(i) = pmb->coarse_x1v(i+1) - pmb->coarse_x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      pmb->coarse_x2v(cjs) = 0.5*(pmb->coarse_x2f(cjs+1) + pmb->coarse_x2f(cjs));
      pmb->coarse_dx2v(cjs) = pmb->coarse_dx2f(cjs);
    } else {
      for (int j=cjs-(pmb->cnghost); j<=cje+(pmb->cnghost); ++j) {
        pmb->coarse_x2v(j) = 0.5*(pmb->coarse_x2f(j+1) + pmb->coarse_x2f(j));
      }
      for (int j=cjs-(pmb->cnghost); j<=cje+(pmb->cnghost)-1; ++j) {
        pmb->coarse_dx2v(j) = pmb->coarse_x2v(j+1) - pmb->coarse_x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      pmb->coarse_x3v(cks) = 0.5*(pmb->coarse_x3f(cks+1) + pmb->coarse_x3f(cks));
      pmb->coarse_dx3v(cks) = pmb->coarse_dx3f(cks);
    } else {
      for (int k=cks-(pmb->cnghost); k<=cke+(pmb->cnghost); ++k) {
        pmb->coarse_x3v(k) = 0.5*(pmb->coarse_x3f(k+1) + pmb->coarse_x3f(k));
      }
      for (int k=cks-(pmb->cnghost); k<=cke+(pmb->cnghost)-1; ++k) {
        pmb->coarse_dx3v(k) = pmb->coarse_x3v(k+1) - pmb->coarse_x3v(k);
      }
    }
  }
}

// destructor

Coordinates::~Coordinates()
{
}

//--------------------------------------------------------------------------------------
// Edge Length functions: returns physical length at cell edges
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))

void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pmy_block->dx1f(i);
  }
  return;
}

// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))

void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pmy_block->dx2f(j);
  }
  return;
}

// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = pmy_block->dx3f(k);
  }
  return;
}

//--------------------------------------------------------------------------------------
// Cell-center Width functions: returns physical width at cell-center

Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return (pmy_block->dx1f(i));
}

Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return (pmy_block->dx2f(j));
}

Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return (pmy_block->dx3f(k));
}

//--------------------------------------------------------------------------------------
// Face Area functions

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area1 = dy dz
    Real& area_i = area(i);
    area_i = (pmy_block->dx2f(j))*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area2 = dx dz
    Real& area_i = area(i);
    area_i = (pmy_block->dx1f(i))*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area3 = dx dy
    Real& area_i = area(i);
    area_i = (pmy_block->dx1f(i))*(pmy_block->dx2f(j));
  }
  return;
}


Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return (pmy_block->dx2f(j))*(pmy_block->dx3f(k));
}

//--------------------------------------------------------------------------------------
// Cell Volume function

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // volume = dx dy dz
    Real& vol_i = vol(i);
    vol_i = (pmy_block->dx1f(i))*(pmy_block->dx2f(j))*(pmy_block->dx3f(k));
  }
  return;
}

Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return pmy_block->dx1f(i)*pmy_block->dx2f(j)*pmy_block->dx3f(k);
}

//--------------------------------------------------------------------------------------
// Coordinate (Geometric) source term functions

void Coordinates::CoordSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}

void Coordinates::CoordSrcTermsX2(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}

void Coordinates::CoordSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  return;
}


//-------------------------------------------------------------------------------------
// Calculate Divv 
void Coordinates::Divv(const AthenaArray<Real> &prim, AthenaArray<Real> &divv)
{

  int is = pmy_block->is; int js = pmy_block->js; int ks = pmy_block->ks;
  int ie = pmy_block->ie; int je = pmy_block->je; int ke = pmy_block->ke;
  int il = is-1; int iu = ie+1;
  int jl, ju, kl, ku;
  Real area_p1, area, vol;
  Real vel_p1, vel;

  if(pmy_block->block_size.nx2 == 1) // 1D
    jl=js, ju=je, kl=ks, ku=ke;
  else if(pmy_block->block_size.nx3 == 1) // 2D
    jl=js-1, ju=je+1, kl=ks, ku=ke;
  else // 3D
    jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;

  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
#pragma simd
      for (int i=il; i<=iu; ++i){
        area_p1 = (pmy_block->dx2f(j))*(pmy_block->dx3f(k));
        area = area_p1; 
        vel_p1 = 0.5*(prim(IM1,k,j,i+1)+prim(IM1,k,j,i));
        vel = 0.5*(prim(IM1,k,j,i)+prim(IM1,k,j,i-1));
        divv(k,j,i) = area_p1*vel_p1 - area*vel;
      }
      if (pmy_block->block_size.nx2 > 1) {
        for (int i=il; i<=iu; ++i){
          area_p1 = (pmy_block->dx1f(i))*(pmy_block->dx3f(k)); 
          area = area_p1;
          vel_p1 = 0.5*(prim(IM2,k,j+1,i)+prim(IM2,k,j,i));
          vel = 0.5*(prim(IM2,k,j,i)+prim(IM2,k,j-1,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      if (pmy_block->block_size.nx3 > 1) {
        for (int i=il; i<=iu; ++i){
          area_p1 = (pmy_block->dx1f(i))*(pmy_block->dx2f(j)); 
          area = area_p1;
          vel_p1 = 0.5*(prim(IM3,k+1,j,i)+prim(IM3,k,j,i));
          vel = 0.5*(prim(IM3,k,j,i)+prim(IM3,k-1,j,i));
          divv(k,j,i) += area_p1*vel_p1 - area*vel;
        }
      }
      for (int i=il; i<=iu; ++i){
        vol = (pmy_block->dx1f(i))*(pmy_block->dx2f(j))*(pmy_block->dx3f(k)); 
        divv(k,j,i) = divv(k,j,i)/vol;
      }
    }
  }

  return;
}


// v_{x1;x1}  covariant derivative at x1 interface
void Coordinates::FaceXdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k,j,i-1))/pmy_block->dx1v(i-1);
  }
  return;
}

// v_{x2;x1}+v_{x1;x2}  covariant derivative at x1 interface
void Coordinates::FaceXdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k,j,i-1))/pmy_block->dx1v(i-1)
             +0.5*(prim(IM1,k,j+1,i)+prim(IM1,k,j+1,i-1)-prim(IM1,k,j-1,i)-prim(IM1,k,j-1,i-1))
              /(pmy_block->dx2v(j-1)+pmy_block->dx2v(j));
  }
  return;
}

// v_{x3;x1}+v_{x1;x3}  covariant derivative at x1 interface
void Coordinates::FaceXdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM3,k,j,i)-prim(IM3,k,j,i-1))/pmy_block->dx1v(i-1)
             +0.5*(prim(IM1,k+1,j,i)+prim(IM1,k+1,j,i-1)-prim(IM1,k-1,j,i)-prim(IM1,k-1,j,i-1))
              /(pmy_block->dx3v(k-1)+pmy_block->dx3v(k));
  }
  return;
}

// v_{x1;x2}+v_{x2;x1}  covariant derivative at x2 interface
void Coordinates::FaceYdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k,j-1,i))/pmy_block->dx2v(j-1)
             +0.5*(prim(IM2,k,j,i+1)+prim(IM2,k,j-1,i+1)-prim(IM2,k,j,i-1)-prim(IM2,k,j-1,i-1))
              /(pmy_block->dx1v(i-1)+pmy_block->dx1v(i));
  }
  return;
}

// v_{x2;x2}  covariant derivative at x2 interface
void Coordinates::FaceYdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k,j-1,i))/pmy_block->dx2v(j-1);
  }
  return;
}

// v_{x3;x2}+v_{x2;x3}  covariant derivative at x2 interface
void Coordinates::FaceYdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM3,k,j,i)-prim(IM3,k,j-1,i))/pmy_block->dx2v(j-1)
             +0.5*(prim(IM2,k+1,j,i)+prim(IM2,k+1,j-1,i)-prim(IM2,k-1,j,i)-prim(IM2,k-1,j-1,i))
              /(pmy_block->dx3v(k-1)+pmy_block->dx3v(k));
  }
  return;
}

// v_{x1;x3}+v_{x3;x1}  covariant derivative at x3 interface
void Coordinates::FaceZdx(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM1,k,j,i)-prim(IM1,k-1,j,i))/pmy_block->dx3v(k-1)
             +0.5*(prim(IM3,k,j,i+1)+prim(IM3,k-1,j,i+1)-prim(IM3,k,j,i-1)-prim(IM3,k-1,j,i-1))
              /(pmy_block->dx1v(i-1)+pmy_block->dx1v(i));
  }
  return;
}

// v_{x2;x3}+v_{x3;x2}  covariant derivative at x3 interface
void Coordinates::FaceZdy(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM2,k,j,i)-prim(IM2,k-1,j,i))/pmy_block->dx3v(k-1)
             +0.5*(prim(IM3,k,j+1,i)+prim(IM3,k-1,j+1,i)-prim(IM3,k,j-1,i)-prim(IM3,k-1,j-1,i))
              /(pmy_block->dx2v(j-1)+pmy_block->dx2v(j));
  } 
  return;
}
    

// v_{x3;x3}  covariant derivative at x3 interface
void Coordinates::FaceZdz(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &prim, AthenaArray<Real> &len)
{   
#pragma simd
  for (int i=il; i<=iu; ++i){
    len(i) = (prim(IM3,k,j,i)-prim(IM3,k-1,j,i))/pmy_block->dx3v(k-1);
  }
  return;
}

//--------------------------------------------------------------------------------------
// Viscous (Geometric) source term functions

void Coordinates::VisSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  return;
}

void Coordinates::VisSrcTermsX2(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  return;
}

void Coordinates::VisSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_p1,
  const AthenaArray<Real> &prim, AthenaArray<Real> &u)
{
  return;
}

