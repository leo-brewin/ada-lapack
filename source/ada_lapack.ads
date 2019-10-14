-- Copyright (c) 1992-2011 The University of Tennessee and The University
--                         of Tennessee Research Foundation.  All rights
--                         reserved.
-- Copyright (c) 2000-2011 The University of California Berkeley. All
--                         rights reserved.
-- Copyright (c) 2006-2012 The University of Colorado Denver.  All rights
--                         reserved.
--
-- $COPYRIGHT$
--
-- Additional copyrights may follow
--
-- $HEADER$
--
-- Redistribution and use in source and binary forms, with or without
-- modification, are permitted provided that the following conditions are
-- met:
--
-- - Redistributions of source code must retain the above copyright
--   notice, this list of conditions and the following disclaimer.
--
-- - Redistributions in binary form must reproduce the above copyright
--   notice, this list of conditions and the following disclaimer listed
--   in this license in the documentation and/or other materials
--   provided with the distribution.
--
-- - Neither the name of the copyright holders nor the names of its
--   contributors may be used to endorse or promote products derived from
--   this software without specific prior written permission.
--
-- The copyright holders provide no reassurances that the source code
-- provided does not infringe any patent, copyright, or any other
-- intellectual property rights of third parties.  The copyright holders
-- disclaim any liability to any recipient for claims brought against
-- recipient by any third party for infringement of that parties
-- intellectual property rights.
--
-- THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
-- "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
-- LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
-- A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
-- OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
-- SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
-- LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
-- DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
-- THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
-- (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
-- OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-- 15-Oct-2012 Version 1.1 Leo Brewin <Leo.Brewin@monash.edu>
-- 15-Aug-2012 Version 1.0 Leo Brewin <Leo.Brewin@monash.edu>

with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Arrays;
with System; use System;

generic

   type Real is digits <>;

   with package Complex_Types  is new Ada.Numerics.Generic_Complex_Types (Real);
   with package Real_Arrays    is new Ada.Numerics.Generic_Real_Arrays (Real);
   with package Complex_Arrays is new Ada.Numerics.Generic_Complex_Arrays (Real_Arrays, Complex_Types);

package Ada_Lapack is

   use Real_Arrays;
   use Complex_Types;
   use Complex_Arrays;

   type Integer_Vector is array (integer range <>) of integer;
      
   lapack_return_code : Integer :=0;

   ----------------------------------------------------------------------------
   -- The publicly visible lapack routines: 
   --
   --   geev, gesv, getrf, getri, getrs, gesdd, and gesvd
   --
   -- The leading letter (d,z) has been dropped from the original name. This is in
   -- part to meet the Lapack licence requirements (that any changes to the Lapack
   -- code must use distinct procedure names) and to consolidate the real and
   -- complex routines under a common (overloaded) procedure name.
   --
   ----------------------------------------------------------------------------

   procedure GEEV
     (JOBVL : Character;
      JOBVR : Character;
      N     : Integer;
      A     : in out Real_Matrix;
      LDA   : Integer;
      WR    : in out Real_Vector;
      WI    : in out Real_Vector;
      VL    : in out Real_Matrix;
      LDVL  : Integer;
      VR    : in out Real_Matrix;
      LDVR  : Integer;
      WORK  : in out Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure GESDD
     (JOBZ  : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Real_Matrix;
      LDA   : Integer;
      S     : in out Real_Vector;
      U     : in out Real_Matrix;
      LDU   : Integer;
      VT    : in out Real_Matrix;
      LDVT  : Integer;
      WORK  : in out Real_Vector;
      LWORK : Integer;
      IWORK : in out Integer_Vector;
      INFO  : in out Integer);

   procedure GESV
     (N    : Integer;
      NRHS : Integer;
      A    : in out Real_Matrix;
      LDA  : Integer;
      IPIV : in out Integer_Vector;
      B    : in out Real_Matrix;
      LDB  : Integer;
      INFO : in out Integer);

   procedure GESVD
     (JOBU  : Character;
      JOBVT : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Real_Matrix;
      LDA   : Integer;
      S     : in out Real_Vector;
      U     : in out Real_Matrix;
      LDU   : Integer;
      VT    : in out Real_Matrix;
      LDVT  : Integer;
      WORK  : in out Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure GETRF
     (M    : Integer;
      N    : Integer;
      A    : in out Real_Matrix;
      LDA  : Integer;
      IPIV : in out Integer_Vector;
      INFO : in out Integer);

   procedure GETRI
     (N     : Integer;
      A     : in out Real_Matrix;
      LDA   : Integer;
      IPIV  : in out Integer_Vector;
      WORK  : in out Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure GETRS
     (TRANS : Character;
      N     : Integer;
      NRHS  : Integer;
      A     : in out Real_Matrix;
      LDA   : Integer;
      IPIV  : in out Integer_Vector;
      B     : in out Real_Matrix;
      LDB   : Integer;
      INFO  : in out Integer);

   procedure SYEV
     (JOBZ  : Character;
      UPLO  : Character;
      N     : Integer;
      A     : in out Real_Matrix;
      LDA   : Integer;
      W     : in out Real_Vector;
      WORK  : in out Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure SYEVD
      (JOBZ   : Character;
       UPLO   : Character;
       N      : Integer;
       A      : in out Real_Matrix;
       LDA    : Integer;
       W      : in out Real_Vector;
       WORK   : in out Real_Vector;
       LWORK  : Integer;
       IWORK  : in out Integer_Vector;
       LIWORK : Integer;
       INFO   : in out Integer);

        procedure SYSV
           (UPLO  : Character;
            N     : Integer;
            NRHS  : Integer;
            A     : in out Real_Matrix;
            LDA   : Integer;
            IPIV  : out    Integer_Vector;
            B     : in out Real_Matrix;
            LDB   : Integer;
            WORK  : in out Real_Vector;
            LWORK : Integer;
            INFO  : in out Integer);

   procedure GEEV
     (JOBVL : Character;
      JOBVR : Character;
      N     : Integer;
      A     : in out Complex_Matrix;
      LDA   : Integer;
      W     : in out Complex_Vector;
      VL    : in out Complex_Matrix;
      LDVL  : Integer;
      VR    : in out Complex_Matrix;
      LDVR  : Integer;
      WORK  : in out Complex_Vector;
      LWORK : Integer;
      RWORK : in out Real_Vector;
      INFO  : in out Integer);

   procedure GESDD
     (JOBZ  : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Complex_Matrix;
      LDA   : Integer;
      S     : in out Real_Vector;
      U     : in out Complex_Matrix;
      LDU   : Integer;
      VT    : in out Complex_Matrix;
      LDVT  : Integer;
      WORK  : in out Complex_Vector;
      LWORK : Integer;
      RWORK : in out Real_Vector;
      IWORK : in out Integer_Vector;
      INFO  : in out Integer);

   procedure GESV
     (N    : Integer;
      NRHS : Integer;
      A    : in out Complex_Matrix;
      LDA  : Integer;
      IPIV : in out Integer_Vector;
      B    : in out Complex_Matrix;
      LDB  : Integer;
      INFO : in out Integer);

   procedure GESVD
     (JOBU  : Character;
      JOBVT : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Complex_Matrix;
      LDA   : Integer;
      S     : in out Real_Vector;
      U     : in out Complex_Matrix;
      LDU   : Integer;
      VT    : in out Complex_Matrix;
      LDVT  : Integer;
      WORK  : in out Complex_Vector;
      LWORK : Integer;
      RWORK : in out Real_Vector;
      INFO  : in out Integer);

   procedure GETRF
     (M    : Integer;
      N    : Integer;
      A    : in out Complex_Matrix;
      LDA  : Integer;
      IPIV : in out Integer_Vector;
      INFO : in out Integer);

   procedure GETRI
     (N     : Integer;
      A     : in out Complex_Matrix;
      LDA   : Integer;
      IPIV  : in out Integer_Vector;
      WORK  : in out Complex_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure GETRS
     (TRANS : Character;
      N     : Integer;
      NRHS  : Integer;
      A     : in out Complex_Matrix;
      LDA   : Integer;
      IPIV  : in out Integer_Vector;
      B     : in out Complex_Matrix;
      LDB   : Integer;
      INFO  : in out Integer);

   procedure HEEV
     (JOBZ  : Character;
      UPLO  : Character;
      N     : Integer;
      A     : in out Complex_Matrix;
      LDA   : Integer;
      W     : in out Real_Vector;
      WORK  : in out Complex_Vector;
      LWORK : Integer;
      RWORK : out Real_Vector;
      INFO  : in out Integer);

   procedure HEEVD
      (JOBZ   : Character;
       UPLO   : Character;
       N      : Integer;
       A      : in out Complex_Matrix;
       LDA    : Integer;
       W      : in out Real_Vector;
       WORK   : in out Complex_Vector;
       LWORK  : Integer;
       RWORK  : in out Real_Vector;
       LRWORK : Integer;
       IWORK  : in out Integer_Vector;
       LIWORK : Integer;
       INFO   : in out Integer)
;

   procedure SYSV
      (UPLO  : Character;
       N     : Integer;
       NRHS  : Integer;
       A     : in out Complex_Matrix;
       LDA   : Integer;
       IPIV  : out    Integer_Vector;
       B     : in out Complex_Matrix;
       LDB   : Integer;
       WORK  : in out Complex_Vector;
       LWORK : Integer;
       INFO  : in out Integer);

private

   type Ftn_Real_Vector is array (Integer range <>) of Real'Base with convention => Fortran;
   type Ftn_Real_Matrix is array (Integer range <>, Integer range <>) of Real'Base with convention => Fortran;
   
   type Ftn_Complex_Vector is array (Integer range <>) of Complex with convention => Fortran;
   type Ftn_Complex_Matrix is array (Integer range <>, Integer range <>) of Complex with convention => Fortran;
   
   type Ftn_Integer_Vector is array (integer range <>) of integer with convention => Fortran;
   type Ftn_Integer_Matrix is array (integer range <>,integer range <>) of integer with convention => Fortran;

   type Ftn_Boolean_Vector is array (integer range <>) of boolean with convention => Fortran;
   type Ftn_Boolean_Matrix is array (integer range <>,integer range <>) of boolean with convention => Fortran;
      
   ----------------------------------------------------------------------------
   -- Blas
   ----------------------------------------------------------------------------

   procedure DAXPY
     (N      : Integer;
      DA     : Real;
      DX_adr : Address;
      INCX   : Integer;
      DY_adr : Address;
      INCY   : Integer);

   procedure DCOPY
     (N      : Integer;
      DX_adr : Address;
      INCX   : Integer;
      DY_adr : Address;
      INCY   : Integer);

   function DDOT
     (N      : Integer;
      DX_adr : Address;
      INCX   : Integer;
      DY_adr : Address;
      INCY   : Integer)
      return   Real;

   procedure DGEMM
     (TRANSA : Character;
      TRANSB : Character;
      M      : Integer;
      N      : Integer;
      K      : Integer;
      ALPHA  : Real;
      A_adr  : Address;
      LDA    : Integer;
      B_adr  : Address;
      LDB    : Integer;
      BETA   : Real;
      C_adr  : Address;
      LDC    : Integer);

   procedure DGEMV
     (TRANS : Character;
      M     : Integer;
      N     : Integer;
      ALPHA : Real;
      A_adr : Address;
      LDA   : Integer;
      X_adr : Address;
      INCX  : Integer;
      BETA  : Real;
      Y_adr : Address;
      INCY  : Integer);

   procedure DGER
     (M     : Integer;
      N     : Integer;
      ALPHA : Real;
      X_adr : Address;
      INCX  : Integer;
      Y_adr : Address;
      INCY  : Integer;
      A_adr : Address;
      LDA   : Integer);

   function DNRM2
     (N     : Integer;
      X_adr : Address;
      INCX  : Integer)
      return  Real;

   procedure DROT
     (N      : Integer;
      DX_adr : Address;
      INCX   : Integer;
      DY_adr : Address;
      INCY   : Integer;
      C      : Real;
      S      : Real);

   procedure DSCAL
     (N      : Integer;
      DA     : Real;
      DX_adr : Address;
      INCX   : Integer);

   procedure DSWAP
     (N      : Integer;
      DX_adr : Address;
      INCX   : Integer;
      DY_adr : Address;
      INCY   : Integer);

   procedure DSYMV
     (UPLO  : Character;
      N     : Integer;
      ALPHA : Real;
      A_adr : Address;
      LDA   : Integer;
      X_adr : Address;
      INCX  : Integer;
      BETA  : Real;
      Y_adr : Address;
      INCY  : Integer);

   procedure DSYR
     (UPLO  : Character;
      N     : Integer;
      ALPHA : Real;
      X_adr : Address;
      INCX  : Integer;
      A_adr : Address;
      LDA   : Integer);

   procedure DSYR2
     (UPLO  : Character;
      N     : Integer;
      ALPHA : Real;
      X_adr : Address;
      INCX  : Integer;
      Y_adr : Address;
      INCY  : Integer;
      A_adr : Address;
      LDA   : Integer);

   procedure DSYR2K
     (UPLO  : Character;
      TRANS : Character;
      N     : Integer;
      K     : Integer;
      ALPHA : Real;
      A_adr : Address;
      LDA   : Integer;
      B_adr : Address;
      LDB   : Integer;
      BETA  : Real;
      C_adr : Address;
      LDC   : Integer);

   procedure DTRMM
     (SIDE   : Character;
      UPLO   : Character;
      TRANSA : Character;
      DIAG   : Character;
      M      : Integer;
      N      : Integer;
      ALPHA  : Real;
      A_adr  : Address;
      LDA    : Integer;
      B_adr  : Address;
      LDB    : Integer);

   procedure DTRMV
     (UPLO  : Character;
      TRANS : Character;
      DIAG  : Character;
      N     : Integer;
      A_adr : Address;
      LDA   : Integer;
      X_adr : Address;
      INCX  : Integer);

   procedure DTRSM
     (SIDE   : Character;
      UPLO   : Character;
      TRANSA : Character;
      DIAG   : Character;
      M      : Integer;
      N      : Integer;
      ALPHA  : Real;
      A_adr  : Address;
      LDA    : Integer;
      B_adr  : Address;
      LDB    : Integer);

   function DZASUM
     (N      : Integer;
      ZX_adr : Address;
      INCX   : Integer)
      return   Real;

   function DZNRM2
     (N     : Integer;
      X_adr : Address;
      INCX  : Integer)
      return  Real;

   function IDAMAX
     (N      : Integer;
      DX_adr : Address;
      INCX   : Integer)
      return   Integer;

   function IZAMAX
     (N      : Integer;
      ZX_adr : Address;
      INCX   : Integer)
      return   Integer;

   procedure ZAXPY
     (N      : Integer;
      ZA     : Complex;
      ZX_adr : Address;
      INCX   : Integer;
      ZY_adr : Address;
      INCY   : Integer);

   procedure ZCOPY
     (N      : Integer;
      ZX_adr : Address;
      INCX   : Integer;
      ZY_adr : Address;
      INCY   : Integer);

   function ZDOTC
     (N      : Integer;
      ZX_adr : Address;
      INCX   : Integer;
      ZY_adr : Address;
      INCY   : Integer)
      return   Complex;

   function ZDOTU
     (N      : Integer;
      ZX_adr : Address;
      INCX   : Integer;
      ZY_adr : Address;
      INCY   : Integer)
      return   Complex;

   procedure ZDROT
     (N      : Integer;
      CX_adr : Address;
      INCX   : Integer;
      CY_adr : Address;
      INCY   : Integer;
      C      : Real;
      S      : Real);

   procedure ZDSCAL
     (N      : Integer;
      DA     : Real;
      ZX_adr : Address;
      INCX   : Integer);

   procedure ZGEMM
     (TRANSA : Character;
      TRANSB : Character;
      M      : Integer;
      N      : Integer;
      K      : Integer;
      ALPHA  : Complex;
      A_adr  : Address;
      LDA    : Integer;
      B_adr  : Address;
      LDB    : Integer;
      BETA   : Complex;
      C_adr  : Address;
      LDC    : Integer);

   procedure ZGEMV
     (TRANS : Character;
      M     : Integer;
      N     : Integer;
      ALPHA : Complex;
      A_adr : Address;
      LDA   : Integer;
      X_adr : Address;
      INCX  : Integer;
      BETA  : Complex;
      Y_adr : Address;
      INCY  : Integer);

   procedure ZGERC
     (M     : Integer;
      N     : Integer;
      ALPHA : Complex;
      X     : in out Ftn_Complex_Vector;
      INCX  : Integer;
      Y     : in out Ftn_Complex_Vector;
      INCY  : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer);

   procedure ZGERU
     (M     : Integer;
      N     : Integer;
      ALPHA : Complex;
      X_adr : Address;
      INCX  : Integer;
      Y_adr : Address;
      INCY  : Integer;
      A_adr : Address;
      LDA   : Integer);

   procedure ZHEMV
     (UPLO  : Character;
      N     : Integer;
      ALPHA : Complex;
      A_adr : Address;
      LDA   : Integer;
      X_adr : Address;
      INCX  : Integer;
      BETA  : Complex;
      Y_adr : Address;
      INCY  : Integer);

   procedure ZHER2
     (UPLO  : Character;
      N     : Integer;
      ALPHA : Complex;
      X_adr : Address;
      INCX  : Integer;
      Y_adr : Address;
      INCY  : Integer;
      A_adr : Address;
      LDA   : Integer);

   procedure ZHER2K
     (UPLO  : Character;
      TRANS : Character;
      N     : Integer;
      K     : Integer;
      ALPHA : Complex;
      A_adr : Address;
      LDA   : Integer;
      B_adr : Address;
      LDB   : Integer;
      BETA  : Real;
      C_adr : Address;
      LDC   : Integer);

   procedure ZSCAL
     (N      : Integer;
      ZA     : Complex;
      ZX_adr : Address;
      INCX   : Integer);

   procedure ZSWAP
     (N      : Integer;
      ZX_adr : Address;
      INCX   : Integer;
      ZY_adr : Address;
      INCY   : Integer);

   procedure ZTRMM
     (SIDE   : Character;
      UPLO   : Character;
      TRANSA : Character;
      DIAG   : Character;
      M      : Integer;
      N      : Integer;
      ALPHA  : Complex;
      A_adr  : Address;
      LDA    : Integer;
      B_adr  : Address;
      LDB    : Integer);

   procedure ZTRMV
     (UPLO  : Character;
      TRANS : Character;
      DIAG  : Character;
      N     : Integer;
      A_adr : Address;
      LDA   : Integer;
      X_adr : Address;
      INCX  : Integer);

   procedure ZTRSM
     (SIDE   : Character;
      UPLO   : Character;
      TRANSA : Character;
      DIAG   : Character;
      M      : Integer;
      N      : Integer;
      ALPHA  : Complex;
      A_adr  : Address;
      LDA    : Integer;
      B_adr  : Address;
      LDB    : Integer);

   procedure ZTRSV
     (UPLO  : Character;
      TRANS : Character;
      DIAG  : Character;
      N     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      X     : in out Ftn_Complex_Vector;
      INCX  : Integer);

   ----------------------------------------------------------------------------
   -- Lapack
   ----------------------------------------------------------------------------

   procedure DBDSDC
     (UPLO      : Character;
      COMPQ     : Character;
      N         : Integer;
      D_adr     : Address;
      E_adr     : Address;
      U_adr     : Address;
      LDU       : Integer;
      VT_adr    : Address;
      LDVT      : Integer;
      Q_adr     : Address;
      IQ_adr    : Address;
      WORK_adr  : Address;
      IWORK_adr : Address;
      INFO      : in out Integer);

   procedure DBDSQR
     (UPLO     : Character;
      N        : Integer;
      NCVT     : Integer;
      NRU      : Integer;
      NCC      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      VT_adr   : Address;
      LDVT     : Integer;
      U_adr    : Address;
      LDU      : Integer;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DGEBAK
     (JOB       : Character;
      SIDE      : Character;
      N         : Integer;
      ILO       : Integer;
      IHI       : Integer;
      SCALE_adr : Address;
      M         : Integer;
      V_adr     : Address;
      LDV       : Integer;
      INFO      : in out Integer);

   procedure DGEBAL
     (JOB       : Character;
      N         : Integer;
      A_adr     : Address;
      LDA       : Integer;
      ILO       : in out Integer;
      IHI       : in out Integer;
      SCALE_adr : Address;
      INFO      : in out Integer);

   procedure DGEBD2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAUQ_adr : Address;
      TAUP_adr : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DGEBRD
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAUQ_adr : Address;
      TAUP_adr : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DGEEV
     (JOBVL : Character;
      JOBVR : Character;
      N     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      WR    : in out Ftn_Real_Vector;
      WI    : in out Ftn_Real_Vector;
      VL    : in out Ftn_Real_Matrix;
      LDVL  : Integer;
      VR    : in out Ftn_Real_Matrix;
      LDVR  : Integer;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DGEHD2
     (N    : Integer;
      ILO  : Integer;
      IHI  : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer;
      TAU  : in out Ftn_Real_Vector;
      WORK : in out Ftn_Real_Vector;
      INFO : in out Integer);

   procedure DGEHRD
     (N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DGELQ2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DGELQF
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DGEQR2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DGEQRF
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DGESDD
     (JOBZ  : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      S     : in out Ftn_Real_Vector;
      U     : in out Ftn_Real_Matrix;
      LDU   : Integer;
      VT    : in out Ftn_Real_Matrix;
      LDVT  : Integer;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      IWORK : in out Ftn_Integer_Vector;
      INFO  : in out Integer);

   procedure DGESV
     (N    : Integer;
      NRHS : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer;
      IPIV : in out Ftn_Integer_Vector;
      B    : in out Ftn_Real_Matrix;
      LDB  : Integer;
      INFO : in out Integer);

   procedure DGESVD
     (JOBU  : Character;
      JOBVT : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      S     : in out Ftn_Real_Vector;
      U     : in out Ftn_Real_Matrix;
      LDU   : Integer;
      VT    : in out Ftn_Real_Matrix;
      LDVT  : Integer;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DGETF2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      IPIV_adr : Address;
      INFO     : in out Integer);

   procedure DGETRF
     (M    : Integer;
      N    : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer;
      IPIV : in out Ftn_Integer_Vector;
      INFO : in out Integer);

   procedure DGETRI
     (N     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      IPIV  : in out Ftn_Integer_Vector;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DGETRS
     (TRANS : Character;
      N     : Integer;
      NRHS  : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      IPIV  : in out Ftn_Integer_Vector;
      B     : in out Ftn_Real_Matrix;
      LDB   : Integer;
      INFO  : in out Integer);

   procedure DHSEQR
     (JOB      : Character;
      COMPZ    : Character;
      N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      H_adr    : Address;
      LDH      : Integer;
      WR_adr   : Address;
      WI_adr   : Address;
      Z_adr    : Address;
      LDZ      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   function DISNAN (DIN : Real) return Boolean;

   procedure DLABAD (SMALL : in out Real; LARGE : in out Real);

   procedure DLABRD
     (M        : Integer;
      N        : Integer;
      NB       : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAUQ_adr : Address;
      TAUP_adr : Address;
      X_adr    : Address;
      LDX      : Integer;
      Y_adr    : Address;
      LDY      : Integer);

   procedure DLACPY
     (UPLO  : Character;
      M     : Integer;
      N     : Integer;
      A_adr : Address;
      LDA   : Integer;
      B_adr : Address;
      LDB   : Integer);

   procedure DLADIV
     (A : Real;
      B : Real;
      C : Real;
      D : Real;
      P : in out Real;
      Q : in out Real);

   procedure DLAE2
     (A   : Real;
      B   : Real;
      C   : Real;
      RT1 : in out Real;
      RT2 : in out Real);

   procedure DLAED0
     (ICOMPQ     : Integer;
      QSIZ       : Integer;
      N          : Integer;
      D_adr      : Address;
      E_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      QSTORE_adr : Address;
      LDQS       : Integer;
      WORK_adr   : Address;
      IWORK_adr  : Address;
      INFO       : in out Integer);

   procedure DLAED1
     (N         : Integer;
      D_adr     : Address;
      Q_adr     : Address;
      LDQ       : Integer;
      INDXQ_adr : Address;
      RHO       : in out Real;
      CUTPNT    : Integer;
      WORK_adr  : Address;
      IWORK_adr : Address;
      INFO      : in out Integer);

   procedure DLAED2
     (K          : in out Integer;
      N          : Integer;
      N1         : Integer;
      D_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      INDXQ_adr  : Address;
      RHO        : in out Real;
      Z_adr      : Address;
      DLAMDA_adr : Address;
      W_adr      : Address;
      Q2_adr     : Address;
      INDX_adr   : Address;
      INDXC_adr  : Address;
      INDXP_adr  : Address;
      COLTYP_adr : Address;
      INFO       : in out Integer);

   procedure DLAED3
     (K          : Integer;
      N          : Integer;
      N1         : Integer;
      D_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      RHO        : Real;
      DLAMDA_adr : Address;
      Q2_adr     : Address;
      INDX_adr   : Address;
      CTOT_adr   : Address;
      W_adr      : Address;
      S_adr      : Address;
      INFO       : in out Integer);

   procedure DLAED4
     (N             : Integer;
      I             : Integer;
      D_adr         : Address;
      Z_adr         : Address;
      Ftn_DELTA_adr : Address;
      RHO           : Real;
      DLAM          : in out Real;
      INFO          : in out Integer);

   procedure DLAED5
     (I             : Integer;
      D_adr         : Address;
      Z_adr         : Address;
      Ftn_DELTA_adr : Address;
      RHO           : Real;
      DLAM          : in out Real);

   procedure DLAED6
     (KNITER : Integer;
      ORGATI : Boolean;
      RHO    : Real;
      D_adr  : Address;
      Z_adr  : Address;
      FINIT  : Real;
      TAU    : in out Real;
      INFO   : in out Integer);

   procedure DLAED7
     (ICOMPQ     : Integer;
      N          : Integer;
      QSIZ       : Integer;
      TLVLS      : Integer;
      CURLVL     : Integer;
      CURPBM     : Integer;
      D_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      INDXQ_adr  : Address;
      RHO        : in out Real;
      CUTPNT     : Integer;
      QSTORE_adr : Address;
      QPTR_adr   : Address;
      PRMPTR_adr : Address;
      PERM_adr   : Address;
      GIVPTR_adr : Address;
      GIVCOL_adr : Address;
      GIVNUM_adr : Address;
      WORK_adr   : Address;
      IWORK_adr  : Address;
      INFO       : in out Integer);

   procedure DLAED8
     (ICOMPQ     : Integer;
      K          : in out Integer;
      N          : Integer;
      QSIZ       : Integer;
      D_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      INDXQ_adr  : Address;
      RHO        : in out Real;
      CUTPNT     : Integer;
      Z_adr      : Address;
      DLAMDA_adr : Address;
      Q2_adr     : Address;
      LDQ2       : Integer;
      W_adr      : Address;
      PERM_adr   : Address;
      GIVPTR     : in out Integer;
      GIVCOL_adr : Address;
      GIVNUM_adr : Address;
      INDXP_adr  : Address;
      INDX_adr   : Address;
      INFO       : in out Integer);

   procedure DLAED9
     (K          : Integer;
      KSTART     : Integer;
      KSTOP      : Integer;
      N          : Integer;
      D_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      RHO        : Real;
      DLAMDA_adr : Address;
      W_adr      : Address;
      S_adr      : Address;
      LDS        : Integer;
      INFO       : in out Integer);

   procedure DLAEDA
     (N          : Integer;
      TLVLS      : Integer;
      CURLVL     : Integer;
      CURPBM     : Integer;
      PRMPTR_adr : Address;
      PERM_adr   : Address;
      GIVPTR_adr : Address;
      GIVCOL_adr : Address;
      GIVNUM_adr : Address;
      Q_adr      : Address;
      QPTR_adr   : Address;
      Z_adr      : Address;
      ZTEMP_adr  : Address;
      INFO       : in out Integer);

   procedure DLAEV2
     (A   : Real;
      B   : Real;
      C   : Real;
      RT1 : in out Real;
      RT2 : in out Real;
      CS1 : in out Real;
      SN1 : in out Real);

   procedure DLAEXC
     (WANTQ : Boolean;
      N     : Integer;
      T     : in out Ftn_Real_Matrix;
      LDT   : Integer;
      Q     : in out Ftn_Real_Matrix;
      LDQ   : Integer;
      J1    : Integer;
      N1    : Integer;
      N2    : Integer;
      WORK  : in out Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure DLAHQR
     (WANTT  : Boolean;
      WANTZ  : Boolean;
      N      : Integer;
      ILO    : Integer;
      IHI    : Integer;
      H_adr  : Address;
      LDH    : Integer;
      WR_adr : Address;
      WI_adr : Address;
      ILOZ   : Integer;
      IHIZ   : Integer;
      Z_adr  : Address;
      LDZ    : Integer;
      INFO   : in out Integer);

   procedure DLAHR2
     (N       : Integer;
      K       : Integer;
      NB      : Integer;
      A_adr   : Address;
      LDA     : Integer;
      TAU_adr : Address;
      T_adr   : Address;
      LDT     : Integer;
      Y_adr   : Address;
      LDY     : Integer);

   function DLAISNAN (DIN1 : Real; DIN2 : Real) return Boolean;

   procedure DLALN2
     (LTRANS : Boolean;
      NA     : Integer;
      NW     : Integer;
      SMIN   : Real;
      CA     : Real;
      A_adr  : Address;
      LDA    : Integer;
      D1     : Real;
      D2     : Real;
      B_adr  : Address;
      LDB    : Integer;
      WR     : Real;
      WI     : Real;
      X_adr  : Address;
      LDX    : Integer;
      SCALE  : in out Real;
      XNORM  : in out Real;
      INFO   : in out Integer);

   function DLAMC3 (A : Real; B : Real) return Real;

   function DLAMCH (CMACH : Character) return Real;

   procedure DLAMRG
     (N1        : Integer;
      N2        : Integer;
      A_adr     : Address;
      DTRD1     : Integer;
      DTRD2     : Integer;
      INDEX_adr : Address);

   function DLANGE
     (NORM : Character;
      M    : Integer;
      N    : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer;
      WORK : in out Ftn_Real_Vector)
      return Real;

   function DLANST
     (NORM  : Character;
      N     : Integer;
      D_adr : Address;
      E_adr : Address)
      return  Real;

   function DLANSY
     (NORM : Character;
      UPLO : Character;
      N    : Integer;
      A    : Ftn_Real_Matrix;
      LDA  : Integer;
      WORK : in out Ftn_Real_Vector)
      return Real;

   procedure DLANV2
     (A    : in out Real;
      B    : in out Real;
      C    : in out Real;
      D    : in out Real;
      RT1R : in out Real;
      RT1I : in out Real;
      RT2R : in out Real;
      RT2I : in out Real;
      CS   : in out Real;
      SN   : in out Real);

   function DLAPY2 (X : Real; Y : Real) return Real;

   function DLAPY3
     (X    : Real;
      Y    : Real;
      Z    : Real)
      return Real;

   procedure DLAQR0
     (WANTT : Boolean;
      WANTZ : Boolean;
      N     : Integer;
      ILO   : Integer;
      IHI   : Integer;
      H     : in out Ftn_Real_Matrix;
      LDH   : Integer;
      WR    : in out Ftn_Real_Vector;
      WI    : in out Ftn_Real_Vector;
      ILOZ  : Integer;
      IHIZ  : Integer;
      Z     : in out Ftn_Real_Matrix;
      LDZ   : Integer;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DLAQR1
     (N     : Integer;
      H_adr : Address;
      LDH   : Integer;
      SR1   : Real;
      SI1   : Real;
      SR2   : Real;
      SI2   : Real;
      V_adr : Address);

   procedure DLAQR2
     (WANTT    : Boolean;
      WANTZ    : Boolean;
      N        : Integer;
      KTOP     : Integer;
      KBOT     : Integer;
      NW       : Integer;
      H_adr    : Address;
      LDH      : Integer;
      ILOZ     : Integer;
      IHIZ     : Integer;
      Z_adr    : Address;
      LDZ      : Integer;
      NS       : in out Integer;
      ND       : in out Integer;
      SR_adr   : Address;
      SI_adr   : Address;
      V_adr    : Address;
      LDV      : Integer;
      NH       : Integer;
      T_adr    : Address;
      LDT      : Integer;
      NV       : Integer;
      WV_adr   : Address;
      LDWV     : Integer;
      WORK_adr : Address;
      LWORK    : Integer);

   procedure DLAQR3
     (WANTT    : Boolean;
      WANTZ    : Boolean;
      N        : Integer;
      KTOP     : Integer;
      KBOT     : Integer;
      NW       : Integer;
      H_adr    : Address;
      LDH      : Integer;
      ILOZ     : Integer;
      IHIZ     : Integer;
      Z_adr    : Address;
      LDZ      : Integer;
      NS       : in out Integer;
      ND       : in out Integer;
      SR_adr   : Address;
      SI_adr   : Address;
      V_adr    : Address;
      LDV      : Integer;
      NH       : Integer;
      T_adr    : Address;
      LDT      : Integer;
      NV       : Integer;
      WV_adr   : Address;
      LDWV     : Integer;
      WORK_adr : Address;
      LWORK    : Integer);

   procedure DLAQR4
     (WANTT    : Boolean;
      WANTZ    : Boolean;
      N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      H_adr    : Address;
      LDH      : Integer;
      WR_adr   : Address;
      WI_adr   : Address;
      ILOZ     : Integer;
      IHIZ     : Integer;
      Z_adr    : Address;
      LDZ      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DLAQR5
     (WANTT  : Boolean;
      WANTZ  : Boolean;
      KACC22 : Integer;
      N      : Integer;
      KTOP   : Integer;
      KBOT   : Integer;
      NSHFTS : Integer;
      SR_adr : Address;
      SI_adr : Address;
      H_adr  : Address;
      LDH    : Integer;
      ILOZ   : Integer;
      IHIZ   : Integer;
      Z_adr  : Address;
      LDZ    : Integer;
      V_adr  : Address;
      LDV    : Integer;
      U_adr  : Address;
      LDU    : Integer;
      NV     : Integer;
      WV_adr : Address;
      LDWV   : Integer;
      NH     : Integer;
      WH_adr : Address;
      LDWH   : Integer);

   procedure DLARF
     (SIDE     : Character;
      M        : Integer;
      N        : Integer;
      V_adr    : Address;
      INCV     : Integer;
      TAU      : Real;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address);

   procedure DLARFB
     (SIDE     : Character;
      TRANS    : Character;
      DIRECT   : Character;
      STOREV   : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      V_adr    : Address;
      LDV      : Integer;
      T_adr    : Address;
      LDT      : Integer;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LDWORK   : Integer);

   procedure DLARFG
     (N     : Integer;
      ALPHA : in out Real;
      X_adr : Address;
      INCX  : Integer;
      TAU   : in out Real);

   procedure DLARFT
     (DIRECT  : Character;
      STOREV  : Character;
      N       : Integer;
      K       : Integer;
      V_adr   : Address;
      LDV     : Integer;
      TAU_adr : Address;
      T_adr   : Address;
      LDT     : Integer);

   procedure DLARFX
     (SIDE     : Character;
      M        : Integer;
      N        : Integer;
      V_adr    : Address;
      TAU      : Real;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address);

   procedure DLARTG
     (F  : Real;
      G  : Real;
      CS : in out Real;
      SN : in out Real;
      R  : in out Real);

   procedure DLAS2
     (F     : Real;
      G     : Real;
      H     : Real;
      SSMIN : in out Real;
      SSMAX : in out Real);

   procedure DLASCL
     (the_TYPE : Character;
      KL       : Integer;
      KU       : Integer;
      CFROM    : Real;
      CTO      : Real;
      M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      INFO     : in out Integer);

   procedure DLASD0
     (N         : Integer;
      SQRE      : Integer;
      D_adr     : Address;
      E_adr     : Address;
      U_adr     : Address;
      LDU       : Integer;
      VT_adr    : Address;
      LDVT      : Integer;
      SMLSIZ    : Integer;
      IWORK_adr : Address;
      WORK_adr  : Address;
      INFO      : in out Integer);

   procedure DLASD1
     (NL        : Integer;
      NR        : Integer;
      SQRE      : Integer;
      D_adr     : Address;
      ALPHA     : in out Real;
      BETA      : in out Real;
      U_adr     : Address;
      LDU       : Integer;
      VT_adr    : Address;
      LDVT      : Integer;
      IDXQ_adr  : Address;
      IWORK_adr : Address;
      WORK_adr  : Address;
      INFO      : in out Integer);

   procedure DLASD2
     (NL         : Integer;
      NR         : Integer;
      SQRE       : Integer;
      K          : in out Integer;
      D_adr      : Address;
      Z_adr      : Address;
      ALPHA      : Real;
      BETA       : Real;
      U_adr      : Address;
      LDU        : Integer;
      VT_adr     : Address;
      LDVT       : Integer;
      DSIGMA_adr : Address;
      U2_adr     : Address;
      LDU2       : Integer;
      VT2_adr    : Address;
      LDVT2      : Integer;
      IDXP_adr   : Address;
      IDX_adr    : Address;
      IDXC_adr   : Address;
      IDXQ_adr   : Address;
      COLTYP_adr : Address;
      INFO       : in out Integer);

   procedure DLASD3
     (NL         : Integer;
      NR         : Integer;
      SQRE       : Integer;
      K          : Integer;
      D_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      DSIGMA_adr : Address;
      U_adr      : Address;
      LDU        : Integer;
      U2_adr     : Address;
      LDU2       : Integer;
      VT_adr     : Address;
      LDVT       : Integer;
      VT2_adr    : Address;
      LDVT2      : Integer;
      IDXC_adr   : Address;
      CTOT_adr   : Address;
      Z_adr      : Address;
      INFO       : in out Integer);

   procedure DLASD4
     (N             : Integer;
      I             : Integer;
      D_adr         : Address;
      Z_adr         : Address;
      Ftn_DELTA_adr : Address;
      RHO           : Real;
      SIGMA         : in out Real;
      WORK_adr      : Address;
      INFO          : in out Integer);

   procedure DLASD5
     (I         : Integer;
      D         : in out Ftn_Real_Vector;
      Z         : in out Ftn_Real_Vector;
      Ftn_DELTA : in out Ftn_Real_Vector;
      RHO       : Real;
      DSIGMA    : in out Real;
      WORK      : in out Ftn_Real_Vector);

   procedure DLASD6
     (ICOMPQ     : Integer;
      NL         : Integer;
      NR         : Integer;
      SQRE       : Integer;
      D_adr      : Address;
      VF_adr     : Address;
      VL_adr     : Address;
      ALPHA      : in out Real;
      BETA       : in out Real;
      IDXQ_adr   : Address;
      PERM_adr   : Address;
      GIVPTR     : in out Integer;
      GIVCOL_adr : Address;
      LDGCOL     : Integer;
      GIVNUM_adr : Address;
      LDGNUM     : Integer;
      POLES_adr  : Address;
      DIFL_adr   : Address;
      DIFR_adr   : Address;
      Z_adr      : Address;
      K          : in out Integer;
      C          : in out Real;
      S          : in out Real;
      WORK_adr   : Address;
      IWORK_adr  : Address;
      INFO       : in out Integer);

   procedure DLASD7
     (ICOMPQ     : Integer;
      NL         : Integer;
      NR         : Integer;
      SQRE       : Integer;
      K          : in out Integer;
      D_adr      : Address;
      Z_adr      : Address;
      ZW_adr     : Address;
      VF_adr     : Address;
      VFW_adr    : Address;
      VL_adr     : Address;
      VLW_adr    : Address;
      ALPHA      : Real;
      BETA       : Real;
      DSIGMA_adr : Address;
      IDX_adr    : Address;
      IDXP_adr   : Address;
      IDXQ_adr   : Address;
      PERM_adr   : Address;
      GIVPTR     : in out Integer;
      GIVCOL_adr : Address;
      LDGCOL     : Integer;
      GIVNUM_adr : Address;
      LDGNUM     : Integer;
      C          : in out Real;
      S          : in out Real;
      INFO       : in out Integer);

   procedure DLASD8
     (ICOMPQ     : Integer;
      K          : Integer;
      D_adr      : Address;
      Z_adr      : Address;
      VF_adr     : Address;
      VL_adr     : Address;
      DIFL_adr   : Address;
      DIFR_adr   : Address;
      LDDIFR     : Integer;
      DSIGMA_adr : Address;
      WORK_adr   : Address;
      INFO       : in out Integer);

   procedure DLASDA
     (ICOMPQ     : Integer;
      SMLSIZ     : Integer;
      N          : Integer;
      SQRE       : Integer;
      D_adr      : Address;
      E_adr      : Address;
      U_adr      : Address;
      LDU        : Integer;
      VT_adr     : Address;
      K_adr      : Address;
      DIFL_adr   : Address;
      DIFR_adr   : Address;
      Z_adr      : Address;
      POLES_adr  : Address;
      GIVPTR_adr : Address;
      GIVCOL_adr : Address;
      LDGCOL     : Integer;
      PERM_adr   : Address;
      GIVNUM_adr : Address;
      C_adr      : Address;
      S_adr      : Address;
      WORK_adr   : Address;
      IWORK_adr  : Address;
      INFO       : in out Integer);

   procedure DLASDQ
     (UPLO     : Character;
      SQRE     : Integer;
      N        : Integer;
      NCVT     : Integer;
      NRU      : Integer;
      NCC      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      VT_adr   : Address;
      LDVT     : Integer;
      U_adr    : Address;
      LDU      : Integer;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DLASDT
     (N         : Integer;
      LVL       : in out Integer;
      ND        : in out Integer;
      INODE_adr : Address;
      NDIML_adr : Address;
      NDIMR_adr : Address;
      MSUB      : Integer);

   procedure DLASET
     (UPLO  : Character;
      M     : Integer;
      N     : Integer;
      ALPHA : Real;
      BETA  : Real;
      A_adr : Address;
      LDA   : Integer);

   procedure DLASQ1
     (N    : Integer;
      D    : in out Ftn_Real_Vector;
      E    : in out Ftn_Real_Vector;
      WORK : in out Ftn_Real_Vector;
      INFO : in out Integer);

   procedure DLASQ2
     (N    : Integer;
      Z    : in out Ftn_Real_Vector;
      INFO : in out Integer);

   procedure DLASQ3
     (I0    : Integer;
      N0    : in out Integer;
      Z     : in out Ftn_Real_Vector;
      PP    : in out Integer;
      DMIN  : in out Real;
      SIGMA : in out Real;
      DESIG : in out Real;
      QMAX  : in out Real; -- doc. says [in] but this routine does modify QMAX, see line 297
      NFAIL : in out Integer;
      ITER  : in out Integer;
      NDIV  : in out Integer;
      IEEE  : Boolean;
      TTYPE : in out Integer;
      DMIN1 : in out Real;
      DMIN2 : in out Real;
      DN    : in out Real;
      DN1   : in out Real;
      DN2   : in out Real;
      G     : in out Real;
      TAU   : in out Real);

   procedure DLASQ4
     (I0    : Integer;
      N0    : Integer;
      Z     : in out Ftn_Real_Vector;
      PP    : Integer;
      N0IN  : Integer;
      DMIN  : Real;
      DMIN1 : Real;
      DMIN2 : Real;
      DN    : Real;
      DN1   : Real;
      DN2   : Real;
      TAU   : in out Real;
      TTYPE : in out Integer;
      G     : in out Real);

   procedure DLASQ5
     (I0    : Integer;
      N0    : Integer;
      Z     : in out Ftn_Real_Vector;
      PP    : Integer;
      TAU   : in out Real;
      SIGMA : Real;
      DMIN  : in out Real;
      DMIN1 : in out Real;
      DMIN2 : in out Real;
      DN    : in out Real;
      DNM1  : in out Real;
      DNM2  : in out Real;
      IEEE  : Boolean;
      EPS   : Real);

   procedure DLASQ6
     (I0    : Integer;
      N0    : Integer;
      Z     : in out Ftn_Real_Vector;
      PP    : Integer;
      DMIN  : in out Real;
      DMIN1 : in out Real;
      DMIN2 : in out Real;
      DN    : in out Real;
      DNM1  : in out Real;
      DNM2  : in out Real);

   procedure DLASR
     (SIDE   : Character;
      PIVOT  : Character;
      DIRECT : Character;
      M      : Integer;
      N      : Integer;
      C_adr  : Address;
      S_adr  : Address;
      A_adr  : Address;
      LDA    : Integer);

   procedure DLASRT
     (ID   : Character;
      N    : Integer;
      D    : in out Ftn_Real_Vector;
      INFO : in out Integer);

   procedure DLASSQ
     (N     : Integer;
      X_adr : Address;
      INCX  : Integer;
      SCALE : in out Real;
      SUMSQ : in out Real);

   procedure DLASV2
     (F     : Real;
      G     : Real;
      H     : Real;
      SSMIN : in out Real;
      SSMAX : in out Real;
      SNR   : in out Real;
      CSR   : in out Real;
      SNL   : in out Real;
      CSL   : in out Real);

   procedure DLASWP
     (N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      K1       : Integer;
      K2       : Integer;
      IPIV_adr : Address;
      INCX     : Integer);

   procedure DLASY2
     (LTRANL : Boolean;
      LTRANR : Boolean;
      ISGN   : Integer;
      N1     : Integer;
      N2     : Integer;
      TL_adr : Address;
      LDTL   : Integer;
      TR_adr : Address;
      LDTR   : Integer;
      B_adr  : Address;
      LDB    : Integer;
      SCALE  : in out Real;
      X_adr  : Address;
      LDX    : Integer;
      XNORM  : in out Real;
      INFO   : in out Integer);

   procedure DLASYF
     (UPLO     : Character;
      N        : Integer;
      NB       : Integer;
      KB       : in out Integer;
      A_adr    : Address;
      LDA      : Integer;
      IPIV_adr : Address;
      W_adr    : Address;
      LDW      : Integer;
      INFO     : in out Integer);

   procedure DLATRD
     (UPLO    : Character;
      N       : Integer;
      NB      : Integer;
      A_adr   : Address;
      LDA     : Integer;
      E_adr   : Address;
      TAU_adr : Address;
      W_adr   : Address;
      LDW     : Integer);

   procedure DORG2L
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DORG2R
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DORGBR
     (VECT     : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORGHR
     (N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORGL2
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DORGLQ
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORGQL
     (M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      TAU   : Ftn_Real_Vector;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DORGQR
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORGTR
     (UPLO     : Character;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORM2L
     (SIDE  : Character;
      TRANS : Character;
      M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      TAU   : Ftn_Real_Vector;
      C     : Ftn_Real_Matrix;
      LDC   : Integer;
      WORK  : Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure DORM2R
     (SIDE  : Character;
      TRANS : Character;
      M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      TAU   : in out Ftn_Real_Vector;
      C     : in out Ftn_Real_Matrix;
      LDC   : Integer;
      WORK  : in out Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure DORMBR
     (VECT     : Character;
      SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORMHR
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORML2
     (SIDE  : Character;
      TRANS : Character;
      M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      TAU   : in out Ftn_Real_Vector;
      C     : in out Ftn_Real_Matrix;
      LDC   : Integer;
      WORK  : in out Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure DORMLQ
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORMQL
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORMQR
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DORMTR
     (SIDE     : Character;
      UPLO     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DSTEDC
     (COMPZ     : Character;
      N         : Integer;
      D_adr     : Address;
      E_adr     : Address;
      Z_adr     : Address;
      LDZ       : Integer;
      WORK_adr  : Address;
      LWORK     : Integer;
      IWORK_adr : Address;
      LIWORK    : Integer;
      INFO      : in out Integer);

   procedure DSTEQR
     (COMPZ    : Character;
      N        : Integer;
      D_adr    : Address;
      E_adr    : Address;
      Z_adr    : Address;
      LDZ      : Integer;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure DSTERF
     (N     : Integer;
      D_adr : Address;
      E_adr : Address;
      INFO  : in out Integer);

   procedure DSYCONV
     (UPLO : Character;
      WAY  : Character;
      N    : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer;
      IPIV : Ftn_Integer_Vector;
      WORK : in out Ftn_Real_Vector;
      INFO : in out Integer);

   procedure DSYEV
     (JOBZ  : Character;
      UPLO  : Character;
      N     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      W     : in out Ftn_Real_Vector;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DSYEVD
     (JOBZ   : Character;
      UPLO   : Character;
      N      : Integer;
      A      : in out Ftn_Real_Matrix;
      LDA    : Integer;
      W      : in out Ftn_Real_Vector;
      WORK   : in out Ftn_Real_Vector;
      LWORK  : Integer;
      IWORK  : in out Ftn_Integer_Vector;
      LIWORK : Integer;
      INFO   : in out Integer);

   procedure DSYSV
     (UPLO  : Character;
      N     : Integer;
      NRHS  : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      IPIV  : out Ftn_Integer_Vector;
      B     : in out Ftn_Real_Matrix;
      LDB   : Integer;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DSYTD2
     (UPLO    : Character;
      N       : Integer;
      A_adr   : Address;
      LDA     : Integer;
      D_adr   : Address;
      E_adr   : Address;
      TAU_adr : Address;
      INFO    : in out Integer);

   procedure DSYTF2
     (UPLO     : Character;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      IPIV_adr : Address;
      INFO     : in out Integer);

   procedure DSYTRD
     (UPLO     : Character;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure DSYTRF
     (UPLO  : Character;
      N     : Integer;
      A     : in out Ftn_Real_Matrix;
      LDA   : Integer;
      IPIV  : in out Ftn_Integer_Vector;
      WORK  : in out Ftn_Real_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure DSYTRS
     (UPLO : Character;
      N    : Integer;
      NRHS : Integer;
      A    : Ftn_Real_Matrix;
      LDA  : Integer;
      IPIV : Ftn_Integer_Vector;
      B    : in out Ftn_Real_Matrix;
      LDB  : Integer;
      INFO : in out Integer);

   procedure DSYTRS2
     (UPLO : Character;
      N    : Integer;
      NRHS : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer;
      IPIV : Ftn_Integer_Vector;
      B    : in out Ftn_Real_Matrix;
      LDB  : Integer;
      WORK : out Ftn_Real_Vector;
      INFO : in out Integer);

   procedure DTREVC
     (SIDE           : Character;
      HOWMNY         : Character;
      the_SELECT_adr : Address;
      N              : Integer;
      T_adr          : Address;
      LDT            : Integer;
      VL_adr         : Address;
      LDVL           : Integer;
      VR_adr         : Address;
      LDVR           : Integer;
      MM             : Integer;
      M              : in out Integer;
      WORK_adr       : Address;
      INFO           : in out Integer);

   procedure DTREXC
     (COMPQ : Character;
      N     : Integer;
      T     : in out Ftn_Real_Matrix;
      LDT   : Integer;
      Q     : in out Ftn_Real_Matrix;
      LDQ   : Integer;
      IFST  : in out Integer;
      ILST  : in out Integer;
      WORK  : in out Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure DTRTI2
     (UPLO  : Character;
      DIAG  : Character;
      N     : Integer;
      A_adr : Address;
      LDA   : Integer;
      INFO  : in out Integer);

   procedure DTRTRI
     (UPLO : Character;
      DIAG : Character;
      N    : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer;
      INFO : in out Integer);

   function IEEECK
     (ISPEC : Integer;
      ZERO  : Real;
      ONE   : Real)
      return  Integer;

   function ILADLC
     (M    : Integer;
      N    : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer)
      return Integer;

   function ILADLR
     (M    : Integer;
      N    : Integer;
      A    : in out Ftn_Real_Matrix;
      LDA  : Integer)
      return Integer;

   function ILAENV
     (ISPEC : Integer;
      NAME  : String;
      OPTS  : String;
      N1    : Integer;
      N2    : Integer;
      N3    : Integer;
      N4    : Integer)
      return  Integer;

   function ILAZLC
     (M    : Integer;
      N    : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer)
      return Integer;

   function ILAZLR
     (M    : Integer;
      N    : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer)
      return Integer;

   function IPARMQ
     (ISPEC : Integer;
      NAME  : String;
      OPTS  : String;
      N     : Integer;
      ILO   : Integer;
      IHI   : Integer;
      LWORK : Integer)
      return  Integer;

   function LSAME (CA : Character; CB : Character) return Boolean;

   procedure XERBLA (SRNAME : String; INFO : Integer);

   procedure ZBDSQR
     (UPLO      : Character;
      N         : Integer;
      NCVT      : Integer;
      NRU       : Integer;
      NCC       : Integer;
      D_adr     : Address;
      E_adr     : Address;
      VT_adr    : Address;
      LDVT      : Integer;
      U_adr     : Address;
      LDU       : Integer;
      C_adr     : Address;
      LDC       : Integer;
      RWORK_adr : Address;
      INFO      : in out Integer);

   procedure ZGEBAK
     (JOB       : Character;
      SIDE      : Character;
      N         : Integer;
      ILO       : Integer;
      IHI       : Integer;
      SCALE_adr : Address;
      M         : Integer;
      V_adr     : Address;
      LDV       : Integer;
      INFO      : in out Integer);

   procedure ZGEBAL
     (JOB       : Character;
      N         : Integer;
      A_adr     : Address;
      LDA       : Integer;
      ILO       : in out Integer;
      IHI       : in out Integer;
      SCALE_adr : Address;
      INFO      : in out Integer);

   procedure ZGEBD2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAUQ_adr : Address;
      TAUP_adr : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure ZGEBRD
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAUQ_adr : Address;
      TAUP_adr : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZGEEV
     (JOBVL : Character;
      JOBVR : Character;
      N     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      W     : in out Ftn_Complex_Vector;
      VL    : in out Ftn_Complex_Matrix;
      LDVL  : Integer;
      VR    : in out Ftn_Complex_Matrix;
      LDVR  : Integer;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      RWORK : in out Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure ZGEHD2
     (N    : Integer;
      ILO  : Integer;
      IHI  : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer;
      TAU  : in out Ftn_Complex_Vector;
      WORK : in out Ftn_Complex_Vector;
      INFO : in out Integer);

   procedure ZGEHRD
     (N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZGELQ2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure ZGELQF
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZGEQR2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure ZGEQRF
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZGESDD
     (JOBZ  : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      S     : in out Ftn_Real_Vector;
      U     : in out Ftn_Complex_Matrix;
      LDU   : Integer;
      VT    : in out Ftn_Complex_Matrix;
      LDVT  : Integer;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      RWORK : in out Ftn_Real_Vector;
      IWORK : in out Ftn_Integer_Vector;
      INFO  : in out Integer);

   procedure ZGESV
     (N    : Integer;
      NRHS : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer;
      IPIV : in out Ftn_Integer_Vector;
      B    : in out Ftn_Complex_Matrix;
      LDB  : Integer;
      INFO : in out Integer);

   procedure ZGESVD
     (JOBU  : Character;
      JOBVT : Character;
      M     : Integer;
      N     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      S     : in out Ftn_Real_Vector;
      U     : in out Ftn_Complex_Matrix;
      LDU   : Integer;
      VT    : in out Ftn_Complex_Matrix;
      LDVT  : Integer;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      RWORK : in out Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure ZGETF2
     (M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      IPIV_adr : Address;
      INFO     : in out Integer);

   procedure ZGETRF
     (M    : Integer;
      N    : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer;
      IPIV : in out Ftn_Integer_Vector;
      INFO : in out Integer);

   procedure ZGETRI
     (N     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      IPIV  : in out Ftn_Integer_Vector;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure ZGETRS
     (TRANS : Character;
      N     : Integer;
      NRHS  : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      IPIV  : in out Ftn_Integer_Vector;
      B     : in out Ftn_Complex_Matrix;
      LDB   : Integer;
      INFO  : in out Integer);

   procedure ZHEEV
     (JOBZ  : Character;
      UPLO  : Character;
      N     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      W     : in out Ftn_Real_Vector;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      RWORK : out Ftn_Real_Vector;
      INFO  : in out Integer);

   procedure ZHEEVD
     (JOBZ   : Character;
      UPLO   : Character;
      N      : Integer;
      A      : in out Ftn_Complex_Matrix;
      LDA    : Integer;
      W      : in out Ftn_Real_Vector;
      WORK   : in out Ftn_Complex_Vector;
      LWORK  : Integer;
      RWORK  : in out Ftn_Real_Vector;
      LRWORK : Integer;
      IWORK  : in out Ftn_Integer_Vector;
      LIWORK : Integer;
      INFO   : in out Integer);

   procedure ZHETD2
     (UPLO    : Character;
      N       : Integer;
      A_adr   : Address;
      LDA     : Integer;
      D_adr   : Address;
      E_adr   : Address;
      TAU_adr : Address;
      INFO    : in out Integer);

   procedure ZHETRD
     (UPLO     : Character;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZHSEQR
     (JOB      : Character;
      COMPZ    : Character;
      N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      H_adr    : Address;
      LDH      : Integer;
      W_adr    : Address;
      Z_adr    : Address;
      LDZ      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZLABRD
     (M        : Integer;
      N        : Integer;
      NB       : Integer;
      A_adr    : Address;
      LDA      : Integer;
      D_adr    : Address;
      E_adr    : Address;
      TAUQ_adr : Address;
      TAUP_adr : Address;
      X_adr    : Address;
      LDX      : Integer;
      Y_adr    : Address;
      LDY      : Integer);

   procedure ZLACGV (N : Integer; X_adr : Address; INCX : Integer);

   procedure ZLACP2
     (UPLO  : Character;
      M     : Integer;
      N     : Integer;
      A_adr : Address;
      LDA   : Integer;
      B_adr : Address;
      LDB   : Integer);

   procedure ZLACPY
     (UPLO  : Character;
      M     : Integer;
      N     : Integer;
      A_adr : Address;
      LDA   : Integer;
      B_adr : Address;
      LDB   : Integer);

   procedure ZLACRM
     (M         : Integer;
      N         : Integer;
      A_adr     : Address;
      LDA       : Integer;
      B_adr     : Address;
      LDB       : Integer;
      C_adr     : Address;
      LDC       : Integer;
      RWORK_adr : Address);

   function ZLADIV (X : Complex; Y : Complex) return Complex;

   procedure ZLAED0
     (QSIZ       : Integer;
      N          : Integer;
      D_adr      : Address;
      E_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      QSTORE_adr : Address;
      LDQS       : Integer;
      RWORK_adr  : Address;
      IWORK_adr  : Address;
      INFO       : in out Integer);

   procedure ZLAED7
     (N          : Integer;
      CUTPNT     : Integer;
      QSIZ       : Integer;
      TLVLS      : Integer;
      CURLVL     : Integer;
      CURPBM     : Integer;
      D_adr      : Address;
      Q_adr      : Address;
      LDQ        : Integer;
      RHO        : in out Real;
      INDXQ_adr  : Address;
      QSTORE_adr : Address;
      QPTR_adr   : Address;
      PRMPTR_adr : Address;
      PERM_adr   : Address;
      GIVPTR_adr : Address;
      GIVCOL_adr : Address;
      GIVNUM_adr : Address;
      WORK_adr   : Address;
      RWORK_adr  : Address;
      IWORK_adr  : Address;
      INFO       : in out Integer);

   procedure ZLAED8
     (K          : in out Integer;
      N          : Integer;
      QSIZ       : Integer;
      Q_adr      : Address;
      LDQ        : Integer;
      D_adr      : Address;
      RHO        : in out Real;
      CUTPNT     : Integer;
      Z_adr      : Address;
      DLAMDA_adr : Address;
      Q2_adr     : Address;
      LDQ2       : Integer;
      W_adr      : Address;
      INDXP_adr  : Address;
      INDX_adr   : Address;
      INDXQ_adr  : Address;
      PERM_adr   : Address;
      GIVPTR     : in out Integer;
      GIVCOL_adr : Address;
      GIVNUM_adr : Address;
      INFO       : in out Integer);

   procedure ZLAHQR
     (WANTT : Boolean;
      WANTZ : Boolean;
      N     : Integer;
      ILO   : Integer;
      IHI   : Integer;
      H_adr : Address;
      LDH   : Integer;
      W_adr : Address;
      ILOZ  : Integer;
      IHIZ  : Integer;
      Z_adr : Address;
      LDZ   : Integer;
      INFO  : in out Integer);

   procedure ZLAHR2
     (N       : Integer;
      K       : Integer;
      NB      : Integer;
      A_adr   : Address;
      LDA     : Integer;
      TAU_adr : Address;
      T_adr   : Address;
      LDT     : Integer;
      Y_adr   : Address;
      LDY     : Integer);

   function ZLANGE
     (NORM : Character;
      M    : Integer;
      N    : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer;
      WORK : in out Ftn_Real_Vector)
      return Real;

   function ZLANHE
     (NORM : Character;
      UPLO : Character;
      N    : Integer;
      A    : Ftn_Complex_Matrix;
      LDA  : Integer;
      WORK : in out Ftn_Real_Vector)
      return Real;

   procedure ZLAQR0
     (WANTT : Boolean;
      WANTZ : Boolean;
      N     : Integer;
      ILO   : Integer;
      IHI   : Integer;
      H     : in out Ftn_Complex_Matrix;
      LDH   : Integer;
      W     : in out Ftn_Complex_Vector;
      ILOZ  : Integer;
      IHIZ  : Integer;
      Z     : in out Ftn_Complex_Matrix;
      LDZ   : Integer;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure ZLAQR1
     (N     : Integer;
      H_adr : Address;
      LDH   : Integer;
      S1    : Complex;
      S2    : Complex;
      V_adr : Address);

   procedure ZLAQR2
     (WANTT    : Boolean;
      WANTZ    : Boolean;
      N        : Integer;
      KTOP     : Integer;
      KBOT     : Integer;
      NW       : Integer;
      H_adr    : Address;
      LDH      : Integer;
      ILOZ     : Integer;
      IHIZ     : Integer;
      Z_adr    : Address;
      LDZ      : Integer;
      NS       : in out Integer;
      ND       : in out Integer;
      SH_adr   : Address;
      V_adr    : Address;
      LDV      : Integer;
      NH       : Integer;
      T_adr    : Address;
      LDT      : Integer;
      NV       : Integer;
      WV_adr   : Address;
      LDWV     : Integer;
      WORK_adr : Address;
      LWORK    : Integer);

   procedure ZLAQR3
     (WANTT    : Boolean;
      WANTZ    : Boolean;
      N        : Integer;
      KTOP     : Integer;
      KBOT     : Integer;
      NW       : Integer;
      H_adr    : Address;
      LDH      : Integer;
      ILOZ     : Integer;
      IHIZ     : Integer;
      Z_adr    : Address;
      LDZ      : Integer;
      NS       : in out Integer;
      ND       : in out Integer;
      SH_adr   : Address;
      V_adr    : Address;
      LDV      : Integer;
      NH       : Integer;
      T_adr    : Address;
      LDT      : Integer;
      NV       : Integer;
      WV_adr   : Address;
      LDWV     : Integer;
      WORK_adr : Address;
      LWORK    : Integer);

   procedure ZLAQR4
     (WANTT    : Boolean;
      WANTZ    : Boolean;
      N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      H_adr    : Address;
      LDH      : Integer;
      W_adr    : Address;
      ILOZ     : Integer;
      IHIZ     : Integer;
      Z_adr    : Address;
      LDZ      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZLAQR5
     (WANTT  : Boolean;
      WANTZ  : Boolean;
      KACC22 : Integer;
      N      : Integer;
      KTOP   : Integer;
      KBOT   : Integer;
      NSHFTS : Integer;
      S_adr  : Address;
      H_adr  : Address;
      LDH    : Integer;
      ILOZ   : Integer;
      IHIZ   : Integer;
      Z_adr  : Address;
      LDZ    : Integer;
      V_adr  : Address;
      LDV    : Integer;
      U_adr  : Address;
      LDU    : Integer;
      NV     : Integer;
      WV_adr : Address;
      LDWV   : Integer;
      NH     : Integer;
      WH_adr : Address;
      LDWH   : Integer);

   procedure ZLARCM
     (M         : Integer;
      N         : Integer;
      A_adr     : Address;
      LDA       : Integer;
      B_adr     : Address;
      LDB       : Integer;
      C_adr     : Address;
      LDC       : Integer;
      RWORK_adr : Address);

   procedure ZLARF
     (SIDE     : Character;
      M        : Integer;
      N        : Integer;
      V_adr    : Address;
      INCV     : Integer;
      TAU      : Complex;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address);

   procedure ZLARFB
     (SIDE     : Character;
      TRANS    : Character;
      DIRECT   : Character;
      STOREV   : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      V_adr    : Address;
      LDV      : Integer;
      T_adr    : Address;
      LDT      : Integer;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LDWORK   : Integer);

   procedure ZLARFG
     (N     : Integer;
      ALPHA : in out Complex;
      X_adr : Address;
      INCX  : Integer;
      TAU   : in out Complex);

   procedure ZLARFT
     (DIRECT  : Character;
      STOREV  : Character;
      N       : Integer;
      K       : Integer;
      V_adr   : Address;
      LDV     : Integer;
      TAU_adr : Address;
      T_adr   : Address;
      LDT     : Integer);

   procedure ZLARTG
     (F  : Complex;
      G  : Complex;
      CS : in out Real;
      SN : in out Complex;
      R  : in out Complex);

   procedure ZLASCL
     (the_TYPE : Character;
      KL       : Integer;
      KU       : Integer;
      CFROM    : Real;
      CTO      : Real;
      M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      INFO     : in out Integer);

   procedure ZLASET
     (UPLO  : Character;
      M     : Integer;
      N     : Integer;
      ALPHA : Complex;
      BETA  : Complex;
      A_adr : Address;
      LDA   : Integer);

   procedure ZLASR
     (SIDE   : Character;
      PIVOT  : Character;
      DIRECT : Character;
      M      : Integer;
      N      : Integer;
      C_adr  : Address;
      S_adr  : Address;
      A_adr  : Address;
      LDA    : Integer);

   procedure ZLASSQ
     (N     : Integer;
      X_adr : Address;
      INCX  : Integer;
      SCALE : in out Real;
      SUMSQ : in out Real);

   procedure ZLASWP
     (N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      K1       : Integer;
      K2       : Integer;
      IPIV_adr : Address;
      INCX     : Integer);

   procedure ZLASYF
     (UPLO     : Character;
      N        : Integer;
      NB       : Integer;
      KB       : in out Integer;
      A_adr    : Address;
      LDA      : Integer;
      IPIV_adr : Address;
      W_adr    : Address;
      LDW      : Integer;
      INFO     : in out Integer);

   procedure ZLATRD
     (UPLO    : Character;
      N       : Integer;
      NB      : Integer;
      A_adr   : Address;
      LDA     : Integer;
      E_adr   : Address;
      TAU_adr : Address;
      W_adr   : Address;
      LDW     : Integer);

   procedure ZLATRS
     (UPLO      : Character;
      TRANS     : Character;
      DIAG      : Character;
      NORMIN    : Character;
      N         : Integer;
      A_adr     : Address;
      LDA       : Integer;
      X_adr     : Address;
      SCALE     : in out Real;
      CNORM_adr : Address;
      INFO      : in out Integer);

   procedure ZROT
     (N      : Integer;
      CX_adr : Address;
      INCX   : Integer;
      CY_adr : Address;
      INCY   : Integer;
      C      : Real;
      S      : Complex);

   procedure ZSTEDC
     (COMPZ     : Character;
      N         : Integer;
      D_adr     : Address;
      E_adr     : Address;
      Z_adr     : Address;
      LDZ       : Integer;
      WORK_adr  : Address;
      LWORK     : Integer;
      RWORK_adr : Address;
      LRWORK    : Integer;
      IWORK_adr : Address;
      LIWORK    : Integer;
      INFO      : in out Integer);

   procedure ZSTEQR
     (COMPZ    : Character;
      N        : Integer;
      D_adr    : Address;
      E_adr    : Address;
      Z_adr    : Address;
      LDZ      : Integer;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure ZSYCONV
     (UPLO : Character;
      WAY  : Character;
      N    : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer;
      IPIV : Ftn_Integer_Vector;
      WORK : in out Ftn_Complex_Vector;
      INFO : in out Integer);

   procedure ZSYR
     (UPLO  : Character;
      N     : Integer;
      ALPHA : Complex;
      X_adr : Address;
      INCX  : Integer;
      A_adr : Address;
      LDA   : Integer);

   procedure ZSYSV
     (UPLO  : Character;
      N     : Integer;
      NRHS  : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      IPIV  : out Ftn_Integer_Vector;
      B     : in out Ftn_Complex_Matrix;
      LDB   : Integer;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure ZSYTF2
     (UPLO     : Character;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      IPIV_adr : Address;
      INFO     : in out Integer);

   procedure ZSYTRF
     (UPLO  : Character;
      N     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      IPIV  : in out Ftn_Integer_Vector;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure ZSYTRS
     (UPLO : Character;
      N    : Integer;
      NRHS : Integer;
      A    : Ftn_Complex_Matrix;
      LDA  : Integer;
      IPIV : Ftn_Integer_Vector;
      B    : in out Ftn_Complex_Matrix;
      LDB  : Integer;
      INFO : in out Integer);

   procedure ZSYTRS2
     (UPLO : Character;
      N    : Integer;
      NRHS : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer;
      IPIV : Ftn_Integer_Vector;
      B    : in out Ftn_Complex_Matrix;
      LDB  : Integer;
      WORK : out Ftn_Complex_Vector;
      INFO : in out Integer);

   procedure ZTREVC
     (SIDE           : Character;
      HOWMNY         : Character;
      the_SELECT_adr : Address;
      N              : Integer;
      T_adr          : Address;
      LDT            : Integer;
      VL_adr         : Address;
      LDVL           : Integer;
      VR_adr         : Address;
      LDVR           : Integer;
      MM             : Integer;
      M              : in out Integer;
      WORK_adr       : Address;
      RWORK_adr      : Address;
      INFO           : in out Integer);

   procedure ZTREXC
     (COMPQ : Character;
      N     : Integer;
      T     : in out Ftn_Complex_Matrix;
      LDT   : Integer;
      Q     : in out Ftn_Complex_Matrix;
      LDQ   : Integer;
      IFST  : Integer;
      ILST  : Integer;
      INFO  : in out Integer);

   procedure ZTRTI2
     (UPLO  : Character;
      DIAG  : Character;
      N     : Integer;
      A_adr : Address;
      LDA   : Integer;
      INFO  : in out Integer);

   procedure ZTRTRI
     (UPLO : Character;
      DIAG : Character;
      N    : Integer;
      A    : in out Ftn_Complex_Matrix;
      LDA  : Integer;
      INFO : in out Integer);

   procedure ZUNG2L
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure ZUNG2R
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure ZUNGBR
     (VECT     : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNGHR
     (N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNGL2
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      INFO     : in out Integer);

   procedure ZUNGLQ
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNGQL
     (M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      TAU   : Ftn_Complex_Vector;
      WORK  : in out Ftn_Complex_Vector;
      LWORK : Integer;
      INFO  : in out Integer);

   procedure ZUNGQR
     (M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNGTR
     (UPLO     : Character;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNM2L
     (SIDE  : Character;
      TRANS : Character;
      M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      TAU   : Ftn_Complex_Vector;
      C     : Ftn_Complex_Matrix;
      LDC   : Integer;
      WORK  : Ftn_Complex_Vector;
      INFO  : in out Integer);

   procedure ZUNM2R
     (SIDE  : Character;
      TRANS : Character;
      M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      TAU   : in out Ftn_Complex_Vector;
      C     : in out Ftn_Complex_Matrix;
      LDC   : Integer;
      WORK  : in out Ftn_Complex_Vector;
      INFO  : in out Integer);

   procedure ZUNMBR
     (VECT     : Character;
      SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNMHR
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      ILO      : Integer;
      IHI      : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNML2
     (SIDE  : Character;
      TRANS : Character;
      M     : Integer;
      N     : Integer;
      K     : Integer;
      A     : in out Ftn_Complex_Matrix;
      LDA   : Integer;
      TAU   : in out Ftn_Complex_Vector;
      C     : in out Ftn_Complex_Matrix;
      LDC   : Integer;
      WORK  : in out Ftn_Complex_Vector;
      INFO  : in out Integer);

   procedure ZUNMLQ
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNMQL
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNMQR
     (SIDE     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      K        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

   procedure ZUNMTR
     (SIDE     : Character;
      UPLO     : Character;
      TRANS    : Character;
      M        : Integer;
      N        : Integer;
      A_adr    : Address;
      LDA      : Integer;
      TAU_adr  : Address;
      C_adr    : Address;
      LDC      : Integer;
      WORK_adr : Address;
      LWORK    : Integer;
      INFO     : in out Integer);

end Ada_Lapack;
