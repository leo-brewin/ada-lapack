generic

package Ada_Lapack.Extras is

   ----------------------------------------------------------------------------------
   -- Some extra procedures built upon the Lapack routines.

	-- determinants of a square marix ------------------------------------------------

   function MatrixDeterm (Source : Real_Matrix) return Real;
   function MatrixDeterm (Source : Complex_Matrix) return Complex;

   -- eigenvalues of a square matrix ------------------------------------------------

   function Eigenvalues (Source : Real_Matrix) return Complex_vector;
   function Eigenvalues (Source : Complex_Matrix) return Complex_vector;

   function EigenvaluesRealSymm (Source : Real_Matrix) return Real_Vector;
   function EigenvaluesHermSymm (Source : Complex_Matrix) return Real_Vector;

   -- eigenvalues and (right) eigenvectors of a square matrix -----------------------

   procedure Eigensystem (Source       :     Real_Matrix;
                          Eigenvalues  : out Complex_Vector;
                          Eigenvectors : out Complex_Matrix);

   procedure Eigensystem (Source       :     Complex_Matrix;
                          Eigenvalues  : out Complex_Vector;
                          Eigenvectors : out Complex_Matrix);

   procedure EigensystemRealSymm (Source       :     Real_Matrix;
                                  Eigenvalues  : out Real_Vector;
                                  Eigenvectors : out Real_Matrix);

   procedure EigensystemHermSymm (Source       :     Complex_Matrix;
                                  Eigenvalues  : out Real_Vector;
                                  Eigenvectors : out Complex_Matrix);

   -- inverse of a square matrix ----------------------------------------------------

   function MatrixInverse (Source : Real_Matrix) return Real_Matrix;
   function MatrixInverse (Source : Complex_Matrix) return Complex_Matrix;

   -- solutions of linear systems of equations --------------------------------------

   function SolveSystem (Source_mat : Real_Matrix;
                         Source_rhs : Real_Vector) return Real_Vector;

   function SolveSystem (Source_mat : Complex_Matrix;
                         Source_rhs : Complex_Vector) return Complex_Vector;

   function SolveSystem (Source_mat : Real_Matrix;
                         Source_rhs : Real_Matrix) return Real_Matrix;

   function SolveSystem (Source_mat : Complex_Matrix;
                         Source_rhs : Complex_Matrix) return Complex_Matrix;

   function SolveSystemRealSymm (Source_mat : Real_Matrix;
                                 Source_rhs : Real_Matrix) return Real_Matrix;

   function SolveSystemHermSymm (Source_mat : Complex_Matrix;
                                 Source_rhs : Complex_Matrix) return Complex_Matrix;

   procedure SolveSystem
      (Solution   : out Real_Vector;
       Source_mat : Real_Matrix;
       Source_rhs : Real_Vector;
       Size       : Integer);

   procedure SolveSystem
      (Solution   : out Complex_Vector;
       Source_mat : Complex_Matrix;
       Source_rhs : Complex_Vector;
       Size       : Integer);

end Ada_Lapack.Extras;
