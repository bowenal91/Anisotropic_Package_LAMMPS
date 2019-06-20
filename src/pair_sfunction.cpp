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

/* ----------------------------------------------------------------------
   Contributing author: Alec Bowen (UChicago)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_sfunction.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

PairSFunction::PairSFunction(LAMMPS *lmp) : Pair(lmp)
{
  //if (lmp->citeme) lmp->citeme->add(cite_pair_gayberne);

  single_enable = 0;
  writedata = 1;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairSFunction::~PairSFunction()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(form);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(sig000);
    memory->destroy(sigcc2);
    memory->destroy(sig220);
    memory->destroy(sig222);
    memory->destroy(sig224);
    memory->destroy(eps000);
    memory->destroy(epscc2);
    memory->destroy(eps220);
    memory->destroy(eps222);
    memory->destroy(eps224);
    memory->destroy(shape1);
    memory->destroy(shape2);
    memory->destroy(well);
    memory->destroy(cut);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    delete [] lshape;
    delete [] setwell;
  }
}

/* ---------------------------------------------------------------------- */

void PairSFunction::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3],neg_r12[3];
  double a1[3][3],b1[3][3],g1[3][3],a2[3][3],b2[3][3],g2[3][3],temp[3][3];
  double fi[3], fj[3], rhat[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *iquat,*jquat;
  double f0,f1,f2;
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    if (form[itype][itype] == ELLIPSE_ELLIPSE) {

      iquat = bonus[ellipsoid[i]].quat;
      MathExtra::quat_to_mat(iquat,a1);
      fi[0] = a1[0][1];
      fi[1] = a1[1][1];
      fi[2] = a1[2][1];
      
    }
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;


      // r12 = center to center vector

      r12[0] = x[i][0]-x[j][0];
      r12[1] = x[i][1]-x[j][1];
      r12[2] = x[i][2]-x[j][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {

        switch (form[itype][jtype]) {
        case SPHERE_SPHERE:
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj *= -r2inv;
          if (eflag) one_eng =
                       r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
                       offset[itype][jtype];
          fforce[0] = -r12[0]*forcelj;
          fforce[1] = -r12[1]*forcelj;
          fforce[2] = -r12[2]*forcelj;
          ttor[0] = ttor[1] = ttor[2] = 0.0;
          rtor[0] = rtor[1] = rtor[2] = 0.0;
          break;

        case SPHERE_ELLIPSE:
          jquat = bonus[ellipsoid[j]].quat;
          MathExtra::quat_to_mat(jquat,a2);
          fj[0] = a2[0][1];
          fj[1] = a2[1][1];
          fj[2] = a2[2][1];
          neg_r12[0] = -r12[0];
          neg_r12[1] = -r12[1];
          neg_r12[2] = -r12[2];
          one_eng = sfunction_lj(j,i,fj,r12,rsq,fforce,rtor);
          ttor[0] = ttor[1] = ttor[2] = 0.0;
          break;

        case ELLIPSE_SPHERE:
          one_eng = sfunction_lj(i,j,fi,r12,rsq,fforce,ttor);
          rtor[0] = rtor[1] = rtor[2] = 0.0;
          break;

        default:
          jquat = bonus[ellipsoid[j]].quat;
          MathExtra::quat_to_mat(jquat,a2);
          fj[0] = a2[0][1];
          fj[1] = a2[1][1];
          fj[2] = a2[2][1];
           
          one_eng = sfunction_analytic(i,j,fi,fj,r12,rsq,
                                      fforce,ttor,rtor);
          break;
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
        tor[i][1] += ttor[1];
        tor[i][2] += ttor[2];

        if (newton_pair || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
          f[j][1] -= fforce[1];
          f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
          tor[j][1] += rtor[1];
          tor[j][2] += rtor[2];
        }

        if (eflag) evdwl = factor_lj*one_eng;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fforce[0],fforce[1],fforce[2],
                                 r12[0],r12[1],r12[2]);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSFunction::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(form,n+1,n+1,"pair:form");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(sig000,n+1,n+1,"pair:sig000");
  memory->create(sigcc2,n+1,n+1,"pair:sigcc2");
  memory->create(sig220,n+1,n+1,"pair:sig220");
  memory->create(sig222,n+1,n+1,"pair:sig222");
  memory->create(sig224,n+1,n+1,"pair:sig224");
  memory->create(eps000,n+1,n+1,"pair:eps000");
  memory->create(epscc2,n+1,n+1,"pair:epscc2");
  memory->create(eps220,n+1,n+1,"pair:eps220");
  memory->create(eps222,n+1,n+1,"pair:eps222");
  memory->create(eps224,n+1,n+1,"pair:eps224");
  memory->create(shape1,n+1,3,"pair:shape1");
  memory->create(shape2,n+1,3,"pair:shape2");
  memory->create(well,n+1,3,"pair:well");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");

  lshape = new double[n+1];
  setwell = new int[n+1];
  for (int i = 1; i <= n; i++) setwell[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSFunction::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSFunction::coeff(int narg, char **arg)
{
  if (narg < 14 || narg > 15)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double eps_one = force->numeric(FLERR,arg[2]);
  double eps000_one = force->numeric(FLERR,arg[3]);
  double epscc2_one = force->numeric(FLERR,arg[4]);
  double eps220_one = force->numeric(FLERR,arg[5]);
  double eps222_one = force->numeric(FLERR,arg[6]);
  double eps224_one = force->numeric(FLERR,arg[7]);
  double sig_one = force->numeric(FLERR,arg[8]);
  double sig000_one = force->numeric(FLERR,arg[9]);
  double sigcc2_one = force->numeric(FLERR,arg[10]);
  double sig220_one = force->numeric(FLERR,arg[11]);
  double sig222_one = force->numeric(FLERR,arg[12]);
  double sig224_one = force->numeric(FLERR,arg[13]);

  
  double cut_one = cut_global;
  if (narg == 15) cut_one = force->numeric(FLERR,arg[14]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = eps_one;
      eps000[i][j] = eps000_one;
      epscc2[i][j] = epscc2_one;
      eps220[i][j] = eps220_one;
      eps222[i][j] = eps222_one;
      eps224[i][j] = eps224_one;
      sigma[i][j] = sig_one;
      sig000[i][j] = sig000_one;
      sigcc2[i][j] = sigcc2_one;
      sig220[i][j] = sig220_one;
      sig222[i][j] = sig222_one;
      sig224[i][j] = sig224_one;

      cut[i][j] = cut_one;
      
      setwell[i] = 1;
      setwell[j] = 1;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSFunction::init_style()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"Pair sfunction requires atom style ellipsoid");

  neighbor->request(this,instance_me);

  // per-type shape precalculations
  // require that atom shapes are identical within each type
  // if shape = 0 for point particle, set shape = 1 as required by Gay-Berne

  for (int i = 1; i <= atom->ntypes; i++) {
    if (!atom->shape_consistency(i,shape1[i][0],shape1[i][1],shape1[i][2]))
      error->all(FLERR,
                 "Pair sfunction requires atoms with same type have same shape");
    if (shape1[i][0] == 0.0)
      shape1[i][0] = shape1[i][1] = shape1[i][2] = 1.0;
    shape2[i][0] = shape1[i][0]*shape1[i][0];
    shape2[i][1] = shape1[i][1]*shape1[i][1];
    shape2[i][2] = shape1[i][2]*shape1[i][2];
    lshape[i] = (shape1[i][0]*shape1[i][1]+shape1[i][2]*shape1[i][2]) *
      sqrt(shape1[i][0]*shape1[i][1]);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSFunction::init_one(int i, int j)
{
  if (setwell[i] == 0 || setwell[j] == 0)
    error->all(FLERR,"Pair sfunction epsilon and sigma coeffs are not all set");

  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  int ishape = 0;
  if (shape1[i][0] != shape1[i][1] ||
      shape1[i][0] != shape1[i][2] ||
      shape1[i][1] != shape1[i][2]) ishape = 1;
  int jshape = 0;
  if (shape1[j][0] != shape1[j][1] ||
      shape1[j][0] != shape1[j][2] ||
      shape1[j][1] != shape1[j][2]) jshape = 1;
  if (ishape == 0 && jshape == 0)
    form[i][i] = form[j][j] = form[i][j] = form[j][i] = SPHERE_SPHERE;
  else if (ishape == 0) {
    form[i][i] = SPHERE_SPHERE; form[j][j] = ELLIPSE_ELLIPSE;
    form[i][j] = SPHERE_ELLIPSE; form[j][i] = ELLIPSE_SPHERE;
  } else if (jshape == 0) {
    form[j][j] = SPHERE_SPHERE; form[i][i] = ELLIPSE_ELLIPSE;
    form[j][i] = SPHERE_ELLIPSE; form[i][j] = ELLIPSE_SPHERE;
  } else
    form[i][i] = form[j][j] = form[i][j] = form[j][i] = ELLIPSE_ELLIPSE;

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  sig000[j][i] = sig000[i][j];
  sigcc2[j][i] = sigcc2[i][j];
  sig220[j][i] = sig220[i][j];
  sig222[j][i] = sig222[i][j];
  sig224[j][i] = sig224[i][j];
  eps000[j][i] = eps000[i][j];
  epscc2[j][i] = epscc2[i][j];
  eps220[j][i] = eps220[i][j];
  eps222[j][i] = eps222[i][j];
  eps224[j][i] = eps224[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSFunction::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    fwrite(&setwell[i],sizeof(int),1,fp);
    if (setwell[i]) fwrite(&well[i][0],sizeof(double),3,fp);
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&eps000[i][j],sizeof(double),1,fp);
        fwrite(&epscc2[i][j],sizeof(double),1,fp);
        fwrite(&eps220[i][j],sizeof(double),1,fp);
        fwrite(&eps222[i][j],sizeof(double),1,fp);
        fwrite(&eps224[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&sig000[i][j],sizeof(double),1,fp);
        fwrite(&sigcc2[i][j],sizeof(double),1,fp);
        fwrite(&sig220[i][j],sizeof(double),1,fp);
        fwrite(&sig222[i][j],sizeof(double),1,fp);
        fwrite(&sig224[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSFunction::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    if (me == 0) fread(&setwell[i],sizeof(int),1,fp);
    MPI_Bcast(&setwell[i],1,MPI_INT,0,world);
    if (setwell[i]) {
      if (me == 0) fread(&well[i][0],sizeof(double),3,fp);
      MPI_Bcast(&well[i][0],3,MPI_DOUBLE,0,world);
    }
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&eps000[i][j],sizeof(double),1,fp);
          fread(&epscc2[i][j],sizeof(double),1,fp);
          fread(&eps220[i][j],sizeof(double),1,fp);
          fread(&eps222[i][j],sizeof(double),1,fp);
          fread(&eps224[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&sig000[i][j],sizeof(double),1,fp);
          fread(&sigcc2[i][j],sizeof(double),1,fp);
          fread(&sig220[i][j],sizeof(double),1,fp);
          fread(&sig222[i][j],sizeof(double),1,fp);
          fread(&sig224[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&eps000[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epscc2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&eps220[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&eps222[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&eps224[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sig000[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigcc2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sig220[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sig222[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sig224[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSFunction::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSFunction::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairSFunction::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g %g %g %g %g %g\n",i,
            epsilon[i][i],eps000[i][i],epscc2[i][i],eps220[i][i],eps222[i][i],eps224[i][i],
            sigma[i][i],sig000[i][i],sigcc2[i][i],sig220[i][i],sig222[i][i],sig224[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairSFunction::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",i,j,
              epsilon[i][j],eps000[i][j],epscc2[i][j],eps220[i][j],eps222[i][j],eps224[i][j],
              sigma[i][j],sig000[i][j],sigcc2[i][j],sig220[i][j],sig222[i][j],sig224[i][j],
              cut[i][j]);
}

/* ----------------------------------------------------------------------
   compute analytic energy, force (fforce), and torque (ttor & rtor)
   based on rotation matrices a and precomputed matrices b and g
   if newton is off, rtor is not calculated for ghost atoms
------------------------------------------------------------------------- */

double PairSFunction::sfunction_analytic(const int i, const int j, double fi[3], double fj[3],
                           double *r12, const double rsqr, double *fforce,
                           double *ttor, double *rtor)
  {
  double tempv[3], tempv2[3];
  double temp[3][3];
  double f0,f1,f2;
  double s0,s000,scc2,s220,s222,s224;
  double e0,e000,ecc2,e220,e222,e224;
  double S000,S202,S022,S220,S222,S224;
  double eps,sig;

  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;

  double r12hat[3];
  MathExtra::normalize3(r12,r12hat);
  double r = sqrt(rsqr);

  // compute the rotational invariants

  f0 = MathExtra::dot3(fi,fj);
  f1 = MathExtra::dot3(fi,r12hat);
  f2 = MathExtra::dot3(fj,r12hat);

  // Specify S-Function Parameters

  s0 = sigma[type[i]][type[j]];
  s000 = sig000[type[i]][type[j]];
  scc2 = sigcc2[type[i]][type[j]];
  s220 = sig220[type[i]][type[j]];
  s222 = sig222[type[i]][type[j]];
  s224 = sig224[type[i]][type[j]];

  e0 = epsilon[type[i]][type[j]];
  e000 = eps000[type[i]][type[j]];
  ecc2 = epscc2[type[i]][type[j]];
  e220 = eps220[type[i]][type[j]];
  e222 = eps222[type[i]][type[j]];
  e224 = eps224[type[i]][type[j]];

  // compute S-Function values
  
  S000 = 1.0;
  S202 = (3.0*f1*f1-1.0)/(2.0*sqrt(5.0));
  S022 = (3.0*f2*f2-1.0)/(2.0*sqrt(5.0));
  S220 = (3.0*f0*f0-1.0)/(2.0*sqrt(5.0)); 
  S222 = (2.0-3.0*f1*f1-3.0*f2*f2-3.0*f0*f0+9.0*f1*f2*f0)/sqrt(70.0);
  S224 = (1.0+2.0*f0*f0-5.0*f1*f1-5.0*f2*f2-20.0*f0*f1*f2+35.0*f1*f1*f2*f2)/(4.0*sqrt(70.0));

  // compute epsilon and sigma for this configuration

  eps = e0*(e000*S000+ecc2*(S202+S022)+e220*S220+e222*S222+e224*S224);
  sig = s0*(s000*S000+scc2*(S202+S022)+s220*S220+s222*S222+s224*S224);
   
  // energy
  // compute u_r
  
  double varrho = s0/(r-sig+s0);
  double varrho6 = pow(varrho,6.0);
  double varrho12 = varrho6*varrho6;
  double u_r = 4.0*eps*(varrho12-varrho6);

  // compute derivatives for force and torque calculations
  
  double dUdr,dUdf0,dUdf1,dUdf2;
  double eps_df0,eps_df1,eps_df2;
  double sig_df0,sig_df1,sig_df2;
    
  S220 = 3.0*f0/sqrt(5.0);
  S222 = (-6.0*f0+9.0*f1*f2)/sqrt(70.0);
  S224 = (4.0*f0-20.0*f1*f2)/(4.0*sqrt(70.0));
  eps_df0 = e0*(e220*S220+e222*S222+e224*S224);
  sig_df0 = s0*(s220*S220+s222*S222+s224*S224);
  
  S202 = 3.0*f1/sqrt(5.0);
  S222 = (-6.0*f1+9.0*f0*f2)/sqrt(70.0);
  S224 = (-10.0*f1-20.0*f0*f2+70.0*f1*f2*f2)/(4.0*sqrt(70.0));
  eps_df1 = e0*(ecc2*S202+e222*S222+e224*S224);
  sig_df1 = s0*(scc2*S202+s222*S222+s224*S224);
  
  S022 = 3*f2/sqrt(5.0);
  S222 = (-6.0*f2+9.0*f0*f1)/sqrt(70.0);
  S224 = (-10.0*f2-20.0*f0*f1+70.0*f1*f1*f2)/(4.0*sqrt(70.0));
  eps_df2 = e0*(ecc2*S022+e222*S222+e224*S224);
  sig_df2 = s0*(scc2*S022+s222*S222+s224*S224);

  dUdr = 4.0*eps*(-12.0*varrho12*varrho + 6.0*varrho6*varrho)/s0;
  dUdf0 = 4.0*eps_df0*(varrho12-varrho6) + 4.0*eps*sig_df0*(12.0*varrho12*varrho-6.0*varrho6*varrho)/s0; 
  dUdf1 = 4.0*eps_df1*(varrho12-varrho6) + 4.0*eps*sig_df1*(12.0*varrho12*varrho-6.0*varrho6*varrho)/s0; 
  dUdf2 = 4.0*eps_df2*(varrho12-varrho6) + 4.0*eps*sig_df2*(12.0*varrho12*varrho-6.0*varrho6*varrho)/s0; 

  // force
 
  fforce[0] = -dUdr * r12hat[0] - dUdf1 * (fi[0]/r - f1*r12[0]/rsqr) - dUdf2 * (fj[0]/r - f2*r12[0]/rsqr);
  fforce[1] = -dUdr * r12hat[1] - dUdf1 * (fi[1]/r - f1*r12[1]/rsqr) - dUdf2 * (fj[1]/r - f2*r12[1]/rsqr);
  fforce[2] = -dUdr * r12hat[2] - dUdf1 * (fi[2]/r - f1*r12[2]/rsqr) - dUdf2 * (fj[2]/r - f2*r12[2]/rsqr);

  // torque for particle 1 and 2

  double temp1[3],temp2[3];

  temp1[0] = -dUdf0*fj[0] - dUdf1*r12hat[0];
  temp1[1] = -dUdf0*fj[1] - dUdf1*r12hat[1];
  temp1[2] = -dUdf0*fj[2] - dUdf1*r12hat[2];

  MathExtra::cross3(fi,temp1,ttor);

  temp2[0] = -dUdf0*fi[0] - dUdf2*r12hat[0]; 
  temp2[1] = -dUdf0*fi[1] - dUdf2*r12hat[1]; 
  temp2[2] = -dUdf0*fi[2] - dUdf2*r12hat[2];

  MathExtra::cross3(fj,temp2,rtor);
 
  return u_r;
}

/* ----------------------------------------------------------------------
   compute analytic energy, force (fforce), and torque (ttor)
   between ellipsoid and lj particle
------------------------------------------------------------------------- */

double PairSFunction::sfunction_lj(const int i, const int j, double fi[3], 
                                   double *r12,const double rsqr, 
                                   double *fforce, double *ttor)

{

  double tempv[3], tempv2[3];
  double temp[3][3];
  double f0,f1,f2;
  double s0,s000,scc2,s220,s222,s224;
  double e0,e000,ecc2,e220,e222,e224;
  double S000,S202,S022,S220,S222,S224;
  double eps,sig;

  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;

  double r12hat[3];
  MathExtra::normalize3(r12,r12hat);
  double r = sqrt(rsqr);

  // compute the rotational invariants

  f1 = MathExtra::dot3(fi,r12hat);

  // Specify S-Function Parameters

  s0 = sigma[type[i]][type[j]];
  s000 = sig000[type[i]][type[j]];
  scc2 = sigcc2[type[i]][type[j]];

  e0 = epsilon[type[i]][type[j]];
  e000 = eps000[type[i]][type[j]];
  ecc2 = epscc2[type[i]][type[j]];

  // compute S-Function values
  
  S000 = 1.0;
  S202 = (3*f1*f1-1)/(2*sqrt(5));

  // compute epsilon and sigma for this configuration

  eps = e0*(e000*S000+ecc2*S202);
  sig = s0*(s000*S000+scc2*S202);

  // energy
  // compute u_r

  double varrho = s0/(r-sig+s0);
  double varrho6 = pow(varrho,6.0);
  double varrho12 = varrho6*varrho6;
  double u_r = 4.0*eps*(varrho12-varrho6);

  // compute derivatives for force and torque calculations

  double dUdr,dUdf1;
  double eps_df1;
  double sig_df1;
    
  S202 = 3*f1/sqrt(5.0);
  eps_df1 = e0*(ecc2*S202);
  sig_df1 = s0*(scc2*S202);
  
  dUdr = 4.0*eps*(-12.0*varrho12*varrho + 6.0*varrho6*varrho)/s0;
  dUdf1 = 4.0*eps_df1*(varrho12-varrho6) + 4.0*eps*sig_df1*(12.0*varrho12*varrho-6.0*varrho6*varrho)/s0; 

  // force
 
  fforce[0] = -dUdr * r12hat[0] - dUdf1 * (fi[0]/r - f1*r12[0]/rsqr);
  fforce[1] = -dUdr * r12hat[1] - dUdf1 * (fi[1]/r - f1*r12[1]/rsqr);
  fforce[2] = -dUdr * r12hat[2] - dUdf1 * (fi[2]/r - f1*r12[2]/rsqr);

  // torque for particle 1 and 2

  double temp1[3];

  temp1[0] = -dUdf1*r12hat[0];
  temp1[1] = -dUdf1*r12hat[1];
  temp1[2] = -dUdf1*r12hat[2];

  MathExtra::cross3(fi,temp1,ttor);

  return u_r;


}

