! -*- cp2k -*-

    &SCF
      EPS_DIIS 0.1
      EPS_SCF 1.0E-4

 @INCLUDE  blank.inc
@inCLude inc_inc.inc

      IGNORE_CONVERGENCE_FAILURE
      SCF_GUESS atomic
    &END SCF
