/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sfunction,PairSFunction)

#else

#ifndef LMP_PAIR_SFUNCTION_H
#define LMP_PAIR_SFUNCTION_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSFunction : public Pair {
 public:
  PairSFunction(LAMMPS *lmp);
  virtual ~PairSFunction();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);

 protected:
  enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

  double cut_global;
  double **cut;

  double gamma,upsilon,mu;   // Gay-Berne parameters
  double **shape1;           // per-type radii in x, y and z
  double **shape2;           // per-type radii in x, y and z SQUARED
  double *lshape;            // precalculation based on the shape
  double **well;             // well depth scaling along each axis ^ -1.0/mu
  double **epsilon,**sigma;  // epsilon and sigma values for atom-type pairs
  double **sig000, **sigcc2, **sig220, **sig222, **sig224; // sigma parameters for S-function potential
  double **eps000, **epscc2, **eps220, **eps222, **eps224; // epsilon parameters for S-function potential

  int **form;
  double **lj1,**lj2,**lj3,**lj4;
  double **offset;
  int *setwell;
  class AtomVecEllipsoid *avec;

  void allocate();
  double sfunction_analytic(const int i, const int j, double fi[3], double fj[3],
                           double *r12, const double rsqr, double *fforce,
                           double *ttor, double *rtor);
  double sfunction_lj(const int i, const int j, double fi[3], 
                           double *r12, const double rsqr,
                           double *fforce, double *ttor);
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair gayberne requires atom style ellipsoid

Self-explanatory.

E: Pair gayberne requires atoms with same type have same shape

Self-explanatory.

E: Pair gayberne epsilon a,b,c coeffs are not all set

Each atom type involved in pair_style gayberne must
have these 3 coefficients set at least once.

E: Bad matrix inversion in mldivide3

This error should not occur unless the matrix is badly formed.

*/
