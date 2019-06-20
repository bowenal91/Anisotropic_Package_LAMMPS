/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "angle_eldihedral.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "math_extra.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05
#define SMALL     0.001
#define SMALLER   0.00001

/* ---------------------------------------------------------------------- */

AngleEldihedral::AngleEldihedral(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleEldihedral::~AngleEldihedral()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k1);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);    
  }
}

/* ---------------------------------------------------------------------- */

void AngleEldihedral::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type,j;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle;
  double dtheta,tk;
  double rsq1,rsq2,s,a,a11,a12,a22;
  double f1[3],f2[3],f3[3],t1[3],t2[3];
  double u1[3],u2[3],a1[3][3],a2[3][3],r1[3],r2[3],r3[3];
  double *quat1,*quat2;
  double c,c11,c22,c33,c12,c13,c23,d12,d13,d23,hi,lo,p,pd;
  double si,siinv,phi;
  double dUdr1[3],dUdr2[3],dUdr3[3],h1[3],h2[3],dummy1[3],dummy2[3];

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"Angle eldihedral requires atom style ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    // 1st bond

    quat1 = bonus[ellipsoid[i1]].quat;
    MathExtra::quat_to_mat(quat1,a1);
    u1[0] = a1[0][1];
    u1[1] = a1[1][1];
    u1[2] = a1[2][1];

    // 1st bond

    r1[0] = -u1[0];
    r1[1] = -u1[1];
    r1[2] = -u1[2];

    // 2nd bond
    
    r2[0] = x[i2][0] - x[i1][0];
    r2[1] = x[i2][1] - x[i1][1];
    r2[2] = x[i2][2] - x[i1][2];

    // 3rd bond

    r3[0] = x[i3][0] - x[i2][0];
    r3[1] = x[i3][1] - x[i2][1];
    r3[2] = x[i3][2] - x[i2][2];
  
    c11 = MathExtra::dot3(r1,r1);
    c22 = MathExtra::dot3(r2,r2);
    c33 = MathExtra::dot3(r3,r3);
    c12 = MathExtra::dot3(r1,r2);
    c13 = MathExtra::dot3(r1,r3);
    c23 = MathExtra::dot3(r2,r3);

    d12 = c11*c22 - c12*c12;
    d13 = c11*c33 - c13*c13;
    d23 = c22*c33 - c23*c23;

    hi = -(c23*c12-c13*c22);
    lo = sqrt(d23*d12);
 
    c = hi/lo;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // force & energy
    // pd = dp/dc

    phi = acos(c);
    si = sin(phi);
    if (fabs(si) < SMALLER) si = SMALLER;
    siinv = 1.0/si;

    p = k1[type]*(1.0 + c) + k2[type]*(1.0 - cos(2.0*phi)) +
      k3[type]*(1.0 + cos(3.0*phi)) + k4[type]*(1.0 - cos(4.0*phi)) ;
    pd = -k1[type] + 2.0*k2[type]*sin(2.0*phi)*siinv -
      3.0*k3[type]*sin(3.0*phi)*siinv + 4.0*k4[type]*sin(4.0*phi)*siinv;

    if (eflag) eangle = p;
    //printf("%f\n",p);
    for (j=0;j<3;j++) {
      dUdr1[j] = (-(c23*r2[j] - c22*r3[j])*lo - hi*(c22*r1[j] - c12*r2[j])*sqrt(d23/d12))/(lo*lo);
      dUdr3[j] = (-(c12*r2[j] - c22*r1[j])*lo - hi*(c22*r3[j] - c23*r2[j])*sqrt(d12/d23))/(lo*lo);

      h1[j] = 2*c33*r2[j] - 2*c23*r3[j];
      h2[j] = 2*c11*r2[j] - 2*c12*r1[j];

      dUdr2[j] = (-(c12*r3[j]+c23*r1[j] - 2*c13*r2[j])*lo - hi*(h1[j]*d12+h2[j]*d23)/(2*sqrt(d12*d23)))/(lo*lo);

    }
    

    for (j=0;j<3;j++) {
      dummy1[j] = dUdr1[j];
      dummy1[j] *= -pd;
      f1[j] = -dUdr2[j];
      f1[j] *= pd;

      f2[j] = dUdr2[j] - dUdr3[j];
      f2[j] *= pd;
      f3[j] = dUdr3[j];
      f3[j] *= pd;

    }

    MathExtra::cross3(u1,dummy1,t1); 
    
    // apply force to each of 4 atoms
      if (newton_bond || i1 < nlocal) {
        f[i1][0] += f1[0];
        f[i1][1] += f1[1];
        f[i1][2] += f1[2];

        tor[i1][0] += t1[0];
        tor[i1][1] += t1[1];
        tor[i1][2] += t1[2];
      }

      if (newton_bond || i2 < nlocal) {
        f[i2][0] += f2[0];
        f[i2][1] += f2[1];
        f[i2][2] += f2[2];
      }

      if (newton_bond || i3 < nlocal) {
        f[i3][0] += f3[0];
        f[i3][1] += f3[1];
        f[i3][2] += f3[2];
      }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleEldihedral::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k1,n+1,"angle:k1");
  memory->create(k2,n+1,"angle:k2");
  memory->create(k3,n+1,"angle:k3");
  memory->create(k4,n+1,"angle:k4");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleEldihedral::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  double k1_one = force->numeric(FLERR,arg[1]);
  double k2_one = force->numeric(FLERR,arg[2]);
  double k3_one = force->numeric(FLERR,arg[3]);
  double k4_one = force->numeric(FLERR,arg[4]);

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k1[i] = k1_one;
    k2[i] = k2_one;
    k3[i] = k3_one;
    k4[i] = k4_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleEldihedral::equilibrium_angle(int i)
{
  //return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleEldihedral::write_restart(FILE *fp)
{
  fwrite(&k1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k3[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k4[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleEldihedral::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k1[1],sizeof(double),atom->nangletypes,fp);
    fread(&k2[1],sizeof(double),atom->nangletypes,fp);
    fread(&k3[1],sizeof(double),atom->nangletypes,fp);
    fread(&k4[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&k1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k4[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleEldihedral::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,k1[i],k2[i],k3[i],k4[i]);
}

/* ---------------------------------------------------------------------- */

double AngleEldihedral::single(int type, int i1, int i2, int i3)
{

    /*
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double dtheta = acos(c) - theta0[type];
  double tk = k[type] * dtheta;
  return tk*dtheta;
  */
}
