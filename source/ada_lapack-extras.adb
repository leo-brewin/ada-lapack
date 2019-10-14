package body Ada_Lapack.Extras is

   function MatrixDeterm (Source : Real_Matrix) return Real
   is

      matrix : Real_Matrix := source;

      matrix_rows : Integer := Matrix'Length (1);
      matrix_cols : Integer := Matrix'Length (2);

      pivots : Integer_Vector (1..matrix_rows);

      determ : Real;

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "MatrixDeterm: source matrix must be square";
      end if;

      GETRF ( A       => matrix,
              M       => matrix_rows,
              N       => matrix_rows,
              LDA     => matrix_rows,
              IPIV    => pivots,
              INFO    => return_code );

      determ := 1.0e0;
      for i in 1..matrix_rows loop
	      determ := determ * matrix(i,i);
         if pivots(i) /= i then
            determ := -determ;
         end if;
      end loop;

      return determ;

   end MatrixDeterm;

   function MatrixDeterm (Source : Complex_Matrix) return Complex
   is

      matrix : Complex_Matrix := source;

      matrix_rows : Integer := Matrix'Length (1);
      matrix_cols : Integer := Matrix'Length (2);

      pivots : Integer_Vector (1..matrix_rows);

      determ : Complex;

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "MatrixDeterm: source matrix must be square";
      end if;

      GETRF ( A       => matrix,
              M       => matrix_rows,
              N       => matrix_rows,
              LDA     => matrix_rows,
              IPIV    => pivots,
              INFO    => return_code );

      determ := Complex'(1.0e0,0.0e0);
      for i in 1..matrix_rows loop
	      determ := determ * matrix(i,i);
         if pivots(i) /= i then
            determ := -determ;
         end if;
      end loop;

      return determ;

   end MatrixDeterm;

   function Eigenvalues (Source : Real_Matrix) return Complex_Vector is

      matrix : Real_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := matrix_rows;

      dummy_Matrix_rows : Integer := matrix_rows;

      dummy_Matrix : Real_Matrix (1 .. dummy_Matrix_rows, 1 .. dummy_Matrix_rows);

      short_vector : Real_Vector (1 .. 1);

      real_eigenvalues : Real_Vector (1 .. eigenvalues_rows);
      imag_eigenvalues : Real_Vector (1 .. eigenvalues_rows);

      eigenvalues : Complex_Vector (1 .. eigenvalues_rows);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[Eigenvalues] source matrix must be square";
      end if;

      GEEV
        (JOBVL => 'N',
         JOBVR => 'N',
         A     => matrix,
         LDA   => matrix_rows,
         N     => matrix_cols,
         WR    => real_eigenvalues,
         WI    => imag_eigenvalues,
         VL    => dummy_Matrix,
         VR    => dummy_Matrix,
         LDVL  => dummy_Matrix_rows,
         LDVR  => dummy_Matrix_rows,
         WORK  => short_vector,
         LWORK => -1,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1));
         work_vector      : Real_Vector (1 .. work_vector_rows);
      begin

         GEEV
           (JOBVL => 'N',
            JOBVR => 'N',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            WR    => real_eigenvalues,
            WI    => imag_eigenvalues,
            VL    => dummy_Matrix,
            VR    => dummy_Matrix,
            LDVL  => dummy_Matrix_rows,
            LDVR  => dummy_Matrix_rows,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[Eigenvalues] failed to find all eigenvalues";
         end if;

         for i in 1 .. eigenvalues_rows loop
            eigenvalues (i) := Complex'(real_eigenvalues (i), imag_eigenvalues (i));
         end loop;

      end;

      return eigenvalues;

   end Eigenvalues;

   function Eigenvalues (Source : Complex_Matrix) return Complex_Vector is

      matrix : Complex_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := matrix_rows;

      dummy_Matrix_rows : Integer := matrix_rows;

      dummy_Matrix : Complex_Matrix (1 .. dummy_Matrix_rows, 1 .. dummy_Matrix_rows);

      short_vector : Complex_Vector (1 .. 1);

      temp_vector : Real_Vector (1 .. 2 * matrix_rows);

      eigenvalues : Complex_Vector (1 .. eigenvalues_rows);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[Eigenvalues] source matrix must be square";
      end if;

      GEEV
        (JOBVL => 'N',
         JOBVR => 'N',
         A     => matrix,
         LDA   => matrix_rows,
         N     => matrix_cols,
         W     => eigenvalues,
         VL    => dummy_Matrix,
         VR    => dummy_Matrix,
         LDVL  => dummy_Matrix_rows,
         LDVR  => dummy_Matrix_rows,
         WORK  => short_vector,
         LWORK => -1,
         RWORK => temp_vector,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1).Re);
         work_vector      : Complex_Vector (1 .. work_vector_rows);
      begin

         GEEV
           (JOBVL => 'N',
            JOBVR => 'N',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            W     => eigenvalues,
            VL    => dummy_Matrix,
            VR    => dummy_Matrix,
            LDVL  => dummy_Matrix_rows,
            LDVR  => dummy_Matrix_rows,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            RWORK => temp_vector,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[Eigenvalues] failed to find all eigenvalues";
         end if;

      end;

      return eigenvalues;

   end Eigenvalues;

   function EigenvaluesRealSymm (Source : Real_Matrix) return Real_Vector is

      matrix : Real_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := matrix_rows;

      eigenvalues : Real_Vector (1 .. eigenvalues_rows);

      short_vector : Real_Vector (1 .. 1);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[EigenvaluesRealSymm] source matrix must be square";
      end if;

      SYEV
        (JOBZ  => 'N',
         UPLO  => 'L',
         A     => matrix,
         N     => matrix_cols,
         LDA   => matrix_rows,
         W     => eigenvalues,
         WORK  => short_vector,
         LWORK => -1,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1));
         work_vector      : Real_Vector (1 .. work_vector_rows);
      begin

         SYEV
           (JOBZ  => 'N',
            UPLO  => 'L',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            W     => eigenvalues,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[EigenvaluesRealSymm] failed to find all eigenvalues";
         end if;

      end;

      return eigenvalues;

   end EigenvaluesRealSymm;

   function EigenvaluesHermSymm (Source : Complex_Matrix) return Real_Vector is

      matrix : Complex_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := matrix_rows;

      eigenvalues : Real_Vector (1 .. eigenvalues_rows);

      short_vector : Complex_Vector (1 .. 1);

      temp_vector : Real_Vector (1 .. 3 * matrix_rows - 2);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[EigenvaluesHermSymm] source matrix must be square";
      end if;

      HEEV
        (JOBZ  => 'N',
         UPLO  => 'L',
         A     => matrix,
         N     => matrix_cols,
         LDA   => matrix_rows,
         W     => eigenvalues,
         WORK  => short_vector,
         LWORK => -1,
         RWORK => temp_vector,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1).Re);
         work_vector      : Complex_Vector (1 .. work_vector_rows);
      begin

         HEEV
           (JOBZ  => 'N',
            UPLO  => 'L',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            W     => eigenvalues,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            RWORK => temp_vector,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[EigenvaluesHermSymm] failed to find all eigenvalues";
         end if;

      end;

      return eigenvalues;

   end EigenvaluesHermSymm;

   procedure Eigensystem
     (Source       : Real_Matrix;
      Eigenvalues  : out Complex_Vector;
      Eigenvectors : out Complex_Matrix)
   is

      i, j : Integer;

      matrix : Real_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := Eigenvalues'Length;

      eigenvectors_rows : Integer := Eigenvectors'Length (1);
      eigenvectors_cols : Integer := Eigenvectors'Length (2);

      real_eigenvalues : Real_Vector (1 .. eigenvalues_rows);
      imag_eigenvalues : Real_Vector (1 .. eigenvalues_rows);

      left_eigenvectors  : Real_Matrix (1 .. eigenvectors_rows, 1 .. eigenvectors_cols);
      right_eigenvectors : Real_Matrix (1 .. eigenvectors_rows, 1 .. eigenvectors_cols);

      short_vector : Real_Vector (1 .. 1);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[Eigensystem] source matrix must be square";
      end if;

      if eigenvectors_rows /= eigenvectors_cols then
         raise Constraint_Error with "[Eigensystem] matrix of eigenvectors must be square";
      end if;

      if eigenvalues_rows /= matrix_rows then
         raise Constraint_Error with "[Eigensystem] matrix and eigenvectors have different no. of rows";
      end if;

      if eigenvectors_rows /= matrix_rows then
         raise Constraint_Error with "[Eigensystem] matrix and eigenvalues have different no. of rows";
      end if;

      GEEV
        (JOBVL => 'N',
         JOBVR => 'V',
         A     => matrix,
         LDA   => matrix_rows,
         N     => matrix_cols,
         WR    => real_eigenvalues,
         WI    => imag_eigenvalues,
         VL    => left_eigenvectors,
         VR    => right_eigenvectors,
         LDVL  => eigenvectors_rows,
         LDVR  => eigenvectors_rows,
         WORK  => short_vector,
         LWORK => -1,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1));
         work_vector      : Real_Vector (1 .. work_vector_rows);
      begin

         GEEV
           (JOBVL => 'N',
            JOBVR => 'V',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            WR    => real_eigenvalues,
            WI    => imag_eigenvalues,
            VL    => left_eigenvectors,
            VR    => right_eigenvectors,
            LDVL  => eigenvectors_rows,
            LDVR  => eigenvectors_rows,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[Eigensystem] failed to find all eigenvector/value pairs";
         end if;

         for i in 1 .. eigenvalues_rows loop
            Eigenvalues (i) := Complex'(real_eigenvalues (i), imag_eigenvalues (i));
         end loop;

         i := 1;
         while i <= eigenvectors_rows loop
            j := 1;
            while j <= eigenvectors_rows loop
               if imag_eigenvalues (j) /= 0.0e0 then
                  Eigenvectors (i, j)     := Complex'(right_eigenvectors (i, j),  right_eigenvectors (i, j + 1));
                  Eigenvectors (i, j + 1) := Complex'(right_eigenvectors (i, j), -right_eigenvectors (i, j + 1));
                  j := j + 2;
               else
                  Eigenvectors (i, j) := Complex'(right_eigenvectors (i, j), 0.0e0);
                  j := j + 1;
               end if;
            end loop;
            i := i + 1;
         end loop;

      end;

   end Eigensystem;

   procedure Eigensystem
     (Source       : Complex_Matrix;
      Eigenvalues  : out Complex_Vector;
      Eigenvectors : out Complex_Matrix)
   is

      matrix : Complex_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := Eigenvalues'Length;

      eigenvectors_rows : Integer := Eigenvectors'Length (1);
      eigenvectors_cols : Integer := Eigenvectors'Length (2);

      left_eigenvectors  : Complex_Matrix (1 .. eigenvectors_rows, 1 .. eigenvectors_cols);
      right_eigenvectors : Complex_Matrix (1 .. eigenvectors_rows, 1 .. eigenvectors_cols);

      short_vector : Complex_Vector (1 .. 1);

      temp_vector : Real_Vector (1 .. 2 * matrix_rows);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[Eigensystem] source matrix must be square";
      end if;

      if eigenvectors_rows /= eigenvectors_cols then
         raise Constraint_Error with "[Eigensystem] matrix of eigenvectors must be square";
      end if;

      if eigenvalues_rows /= matrix_rows then
         raise Constraint_Error with "[Eigensystem] matrix and eigenvectors have different no. of rows";
      end if;

      if eigenvectors_rows /= matrix_rows then
         raise Constraint_Error with "[Eigensystem] matrix and eigenvalues have different no. of rows";
      end if;

      GEEV
        (JOBVL => 'N',
         JOBVR => 'V',
         A     => matrix,
         LDA   => matrix_rows,
         N     => matrix_cols,
         W     => Eigenvalues,
         VL    => left_eigenvectors,
         VR    => right_eigenvectors,
         LDVL  => eigenvectors_rows,
         LDVR  => eigenvectors_rows,
         WORK  => short_vector,
         LWORK => -1,
         RWORK => temp_vector,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1).Re);
         work_vector      : Complex_Vector (1 .. work_vector_rows);
      begin

         GEEV
           (JOBVL => 'N',
            JOBVR => 'V',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            W     => Eigenvalues,
            VL    => left_eigenvectors,
            VR    => right_eigenvectors,
            LDVL  => eigenvectors_rows,
            LDVR  => eigenvectors_rows,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            RWORK => temp_vector,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[Eigensystem] failed to find all eigenvectors/values pairs";
         end if;

         Eigenvectors := right_eigenvectors;

      end;

   end Eigensystem;

   procedure EigensystemRealSymm
     (Source       : Real_Matrix;
      Eigenvalues  : out Real_Vector;
      Eigenvectors : out Real_Matrix)
   is

      matrix : Real_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := Eigenvalues'Length;

      eigenvectors_rows : Integer := Eigenvectors'Length (1);
      eigenvectors_cols : Integer := Eigenvectors'Length (2);

      short_vector : Real_Vector (1 .. 1);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[EigensystemRealSymm] source matrix must be square";
      end if;

      if eigenvectors_rows /= eigenvectors_cols then
         raise Constraint_Error with "[EigensystemRealSymm] matrix of eigenvectors must be square";
      end if;

      if eigenvalues_rows /= matrix_rows then
         raise Constraint_Error with "[EigensystemRealSymm] matrix and eigenvectors have different no. of rows";
      end if;

      if eigenvectors_rows /= matrix_rows then
         raise Constraint_Error with "[EigensystemRealSymm] matrix and eigenvalues have different no. of rows";
      end if;

      SYEV
        (JOBZ  => 'V',
         UPLO  => 'L',
         A     => matrix,
         N     => matrix_cols,
         LDA   => matrix_rows,
         W     => Eigenvalues,
         WORK  => short_vector,
         LWORK => -1,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1));
         work_vector      : Real_Vector (1 .. work_vector_rows);
      begin

         SYEV
           (JOBZ  => 'V',
            UPLO  => 'L',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            W     => Eigenvalues,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[EigensystemRealSymm] failed to find all eigenvector/value pairs";
         end if;

         for i in 1 .. matrix_rows loop
            for j in 1 .. matrix_cols loop
               Eigenvectors (i, j) := matrix (i, j);
            end loop;
         end loop;

      end;

   end EigensystemRealSymm;

   procedure EigensystemHermSymm
     (Source       : Complex_Matrix;
      Eigenvalues  : out Real_Vector;
      Eigenvectors : out Complex_Matrix)
   is

      matrix : Complex_Matrix := Source;

      matrix_rows : Integer := matrix'Length (1);
      matrix_cols : Integer := matrix'Length (2);

      eigenvalues_rows : Integer := Eigenvalues'Length;

      eigenvectors_rows : Integer := Eigenvectors'Length (1);
      eigenvectors_cols : Integer := Eigenvectors'Length (2);

      short_vector : Complex_Vector (1 .. 1);

      temp_vector : Real_Vector (1 .. 3 * matrix_rows - 2);

      return_code : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "[EigensystemHermSymm] source matrix must be square";
      end if;

      if eigenvectors_rows /= eigenvectors_cols then
         raise Constraint_Error with "[EigensystemHermSymm] matrix of eigenvectors must be square";
      end if;

      if eigenvalues_rows /= matrix_rows then
         raise Constraint_Error with "[EigensystemHermSymm] matrix and eigenvectors have different no. of rows";
      end if;

      if eigenvectors_rows /= matrix_rows then
         raise Constraint_Error with "[EigensystemHermSymm] matrix and eigenvalues have different no. of rows";
      end if;

      HEEV
        (JOBZ  => 'V',
         UPLO  => 'L',
         A     => matrix,
         N     => matrix_cols,
         LDA   => matrix_rows,
         W     => Eigenvalues,
         WORK  => short_vector,
         LWORK => -1,
         RWORK => temp_vector,
         INFO  => return_code);

      declare
         work_vector_rows : Integer := Integer (short_vector (1).Re);
         work_vector      : Complex_Vector (1 .. work_vector_rows);
      begin

         HEEV
           (JOBZ  => 'V',
            UPLO  => 'L',
            A     => matrix,
            N     => matrix_cols,
            LDA   => matrix_rows,
            W     => Eigenvalues,
            WORK  => work_vector,
            LWORK => work_vector_rows,
            RWORK => temp_vector,
            INFO  => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "[EigensystemHermSymm] failed to find all eigenvectors/values pairs";
         end if;

         for i in 1 .. matrix_rows loop
            for j in 1 .. matrix_cols loop
               Eigenvectors (i, j) := matrix (i, j);
            end loop;
         end loop;

      end;

   end EigensystemHermSymm;

   function MatrixInverse (Source : Real_Matrix) return Real_Matrix
   is

      matrix : Real_Matrix := source;

      matrix_rows : Integer := Matrix'Length (1);
      matrix_cols : Integer := Matrix'Length (2);

      pivots : Integer_Vector (1..matrix_rows);

      short_vector : Real_Vector (1..1);

      return_code  : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "MatrixInverse: source matrix must be square";
      end if;

      GETRF ( A       => matrix,
              M       => matrix_rows,
              N       => matrix_rows,
              LDA     => matrix_rows,
              IPIV    => pivots,
              INFO    => return_code );

      GETRI ( A       => matrix,
              N       => matrix_rows,
              LDA     => matrix_rows,
              IPIV    => pivots,
              WORK    => short_vector,
              LWORK   => -1,
              INFO    => return_code);

      declare
         work_vector_rows : Integer := Integer( short_vector(1) );
         work_vector : Real_Vector (1 .. work_vector_rows);
      begin

         GETRI ( A       => matrix,
                 N       => matrix_rows,
                 LDA     => matrix_rows,
                 IPIV    => pivots,
                 WORK    => work_vector,
                 LWORK   => work_vector_rows,
                 INFO    => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "MatrixInverse: matrix is singular";
         end if;

      end;

      return matrix;

   end MatrixInverse;

   function MatrixInverse (Source : Complex_Matrix) return Complex_Matrix
   is

      matrix : Complex_Matrix := source;

      matrix_rows : Integer := Matrix'Length (1);
      matrix_cols : Integer := Matrix'Length (2);

      pivots : Integer_Vector (1..matrix_rows);

      short_vector : Complex_Vector (1..1);

      return_code  : Integer := 0;

   begin

      if matrix_rows /= matrix_cols then
         raise Constraint_Error with "MatrixInverse: source matrix must be square";
      end if;

      GETRF ( A       => matrix,
              M       => matrix_rows,
              N       => matrix_rows,
              LDA     => matrix_rows,
              IPIV    => pivots,
              INFO    => return_code );

      GETRI ( A       => matrix,
              N       => matrix_rows,
              LDA     => matrix_rows,
              IPIV    => pivots,
              WORK    => short_vector,
              LWORK   => -1,
              INFO    => return_code);

      declare
         work_vector_rows : Integer := Integer( short_vector(1).Re );
         work_vector : Complex_Vector (1 .. work_vector_rows);
      begin

         GETRI ( A       => matrix,
                 N       => matrix_rows,
                 LDA     => matrix_rows,
                 IPIV    => pivots,
                 WORK    => work_vector,
                 LWORK   => work_vector_rows,
                 INFO    => return_code);

         if return_code /= 0 then
            lapack_return_code := return_code;
            raise Constraint_Error with "MatrixInverse: matrix is singular";
         end if;

      end;

      return matrix;

   end MatrixInverse;

	function SolveSystem
	   (Source_mat : Real_Matrix;
	    Source_rhs : Real_Vector)
	    return       Real_Vector
	is

	   matrix        : Real_Matrix := Source_mat;
	   rhs           : Real_Vector := Source_rhs;

	   matrix_rows   : Integer := matrix'Length (1);
	   matrix_cols   : Integer := matrix'Length (2);

	   rhs_rows      : Integer := rhs'Length (1);
	   rhs_cols      : Integer := 1;

	   mixed_sol_rhs : Real_Matrix (1..matrix_rows,1..1);
	   solution      : Real_Vector (1..matrix_rows);

	   pivots        : Integer_Vector (1..matrix_rows);

	   return_code   : Integer := 0;

	begin

	   if matrix_rows /= matrix_cols then
	      raise Constraint_Error with "[SolveSystem] source matrix must be square";
	   end if;

	   if rhs_rows /= matrix_rows then
	      raise Constraint_Error with "[SolveSystem] matrix and rhs must have equal no. of rows";
	   end if;

	   for i in 1..rhs_rows loop
	      mixed_sol_rhs (i,1) := rhs (i);
	   end loop;

	   GESV (A    => matrix,
	         LDA  => matrix_rows,
	         N    => matrix_rows,
	         IPIV => pivots,
	         B    => mixed_sol_rhs,
	         LDB  => matrix_rows,
	         NRHS => 1,
	         INFO => return_code);

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystem] matrix may be singular";
	   end if;

	   for i in 1..rhs_rows loop
	      solution (i) := mixed_sol_rhs (i,1);
	   end loop;

	   return solution;

	end SolveSystem;

	function SolveSystem
	  (Source_mat : Complex_Matrix;
	   Source_rhs : Complex_Vector)
	   return       Complex_Vector
	is

	   matrix        : Complex_Matrix := Source_mat;
	   rhs           : Complex_Vector := Source_rhs;

	   matrix_rows   : Integer := matrix'Length (1);
	   matrix_cols   : Integer := matrix'Length (2);

	   rhs_rows      : Integer := rhs'Length (1);
	   rhs_cols      : Integer := 1;

	   mixed_sol_rhs : Complex_Matrix (1..matrix_rows,1..1);
	   solution      : Complex_Vector (1..matrix_rows);

	   pivots        : Integer_Vector (1..matrix_rows);

	   return_code   : Integer := 0;

	begin

	   if matrix_rows /= matrix_cols then
	      raise Constraint_Error with "[SolveSystem] source matrix must be square";
	   end if;

	   if rhs_rows /= matrix_rows then
	      raise Constraint_Error with "[SolveSystem] matrix and rhs must have equal no. of rows";
	   end if;

	   for i in 1..rhs_rows loop
	      mixed_sol_rhs (i,1) := rhs (i);
	   end loop;

	   GESV (A    => matrix,
	         LDA  => matrix_rows,
	         N    => matrix_rows,
	         IPIV => pivots,
	         B    => mixed_sol_rhs,
	         LDB  => matrix_rows,
	         NRHS => 1,
	         INFO => return_code);

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystem] matrix may be singular";
	   end if;

	   for i in 1..rhs_rows loop
	      solution (i) := mixed_sol_rhs (i,1);
	   end loop;

	   return solution;

	end SolveSystem;

	function SolveSystem
	  (Source_mat : Real_Matrix;
	   Source_rhs : Real_Matrix)
	   return       Real_Matrix
	is

	   matrix        : Real_Matrix := Source_mat;
	   rhs           : Real_Matrix := Source_rhs;

	   matrix_rows   : Integer := matrix'Length (1);
	   matrix_cols   : Integer := matrix'Length (2);

	   rhs_rows      : Integer := rhs'Length (1);
	   rhs_cols      : Integer := rhs'Length (2);

	   solution      : Real_Matrix := rhs;

	   solution_rows : Integer := rhs_rows;
	   solution_cols : Integer := rhs_cols;

	   pivots        : Integer_Vector (1 .. matrix_rows);

	   return_code   : Integer := 0;

	begin

	   if matrix_rows /= matrix_cols then
	      raise Constraint_Error with "[SolveSystem] source matrix must be square";
	   end if;

	   if rhs_rows /= matrix_rows then
	      raise Constraint_Error with "[SolveSystem] matrix and rhs must have equal no. of rows";
	   end if;

	   GESV (A    => matrix,
	         LDA  => matrix_rows,
	         N    => matrix_rows,
	         IPIV => pivots,
	         B    => solution,
	         LDB  => solution_rows,
	         NRHS => solution_cols,
	         INFO => return_code);

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystem] matrix may be singular";
	   end if;

	   return solution;

	end SolveSystem;

	function SolveSystem
	  (Source_mat : Complex_Matrix;
	   Source_rhs : Complex_Matrix)
	   return       Complex_Matrix
	is

	   matrix        : Complex_Matrix := Source_mat;
	   rhs           : Complex_Matrix := Source_rhs;

	   matrix_rows   : Integer := matrix'Length (1);
	   matrix_cols   : Integer := matrix'Length (2);

	   rhs_rows      : Integer := rhs'Length (1);
	   rhs_cols      : Integer := rhs'Length (2);

	   solution      : Complex_Matrix := rhs;

	   solution_rows : Integer := rhs_rows;
	   solution_cols : Integer := rhs_cols;

	   pivots        : Integer_Vector (1 .. matrix_rows);

	   return_code   : Integer := 0;

	begin

	   if matrix_rows /= matrix_cols then
	      raise Constraint_Error with "[SolveSystem] source matrix must be square";
	   end if;

	   if rhs_rows /= matrix_rows then
	      raise Constraint_Error with "[SolveSystem] matrix and rhs must have equal no. of rows";
	   end if;

	   GESV (A    => matrix,
	         LDA  => matrix_rows,
	         N    => matrix_rows,
	         IPIV => pivots,
	         B    => solution,
	         LDB  => solution_rows,
	         NRHS => solution_cols,
	         INFO => return_code);

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystem] matrix may be singular";
	   end if;

	   return solution;

	end SolveSystem;

	function SolveSystemRealSymm
	  (Source_mat : Real_Matrix;
	   Source_rhs : Real_Matrix)
	   return       Real_Matrix
	is

	   matrix       : Real_Matrix := Source_mat;
	   rhs          : Real_Matrix := Source_rhs;

	   matrix_rows  : Integer := matrix'Length (1);
	   matrix_cols  : Integer := matrix'Length (2);

	   rhs_rows     : Integer := rhs'Length (1);
	   rhs_cols     : Integer := rhs'Length (2);

	   pivots       : Integer_Vector (1 .. matrix_rows);

	   short_vector : Real_Vector (1..1);

	   return_code  : Integer := 0;

	begin

	   if matrix_rows /= matrix_cols then
	      raise Constraint_Error with "[SolveSystemRealSymm] source matrix must be square";
	   end if;

	   if rhs_rows /= matrix_rows then
	      raise Constraint_Error with "[SolveSystemRealSymm] matrix and rhs must have equal no. of rows";
	   end if;

	   SYSV ( UPLO  => 'L',
	          A     => matrix,
	          LDA   => matrix_rows,
	          N     => matrix_cols,
	          B     => rhs,
	          LDB   => rhs_rows,
	          NRHS  => rhs_cols,
	          IPIV  => pivots,
	          WORK  => short_vector,
	          LWORK => -1,
	          INFO  => return_code );

	   declare
	      work_vector_max : Constant Integer := Integer( short_vector(1) );
	      work_vector     : Real_Vector (1 .. work_vector_max);
	   begin

	      SYSV ( UPLO  => 'L',
	             A     => matrix,
	             LDA   => matrix_rows,
	             N     => matrix_cols,
	             B     => rhs,
	             LDB   => rhs_rows,
	             NRHS  => rhs_cols,
	             IPIV  => pivots,
	             WORK  => work_vector,
	             LWORK => work_vector_max,
	             INFO  => return_code );

	    end;

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystemRealSymm] matrix may be singular";
	   end if;

	   return rhs;

	end SolveSystemRealSymm;

	function SolveSystemHermSymm
	  (Source_mat : Complex_Matrix;
	   Source_rhs : Complex_Matrix)
	   return       Complex_Matrix
	is

	   matrix       : Complex_Matrix := Source_mat;
	   rhs          : Complex_Matrix := Source_rhs;

	   matrix_rows  : Integer := matrix'Length (1);
	   matrix_cols  : Integer := matrix'Length (2);

	   rhs_rows     : Integer := rhs'Length (1);
	   rhs_cols     : Integer := rhs'Length (2);

	   pivots       : Integer_Vector (1 .. matrix_rows);

	   short_vector : Complex_Vector (1..1);

	   return_code  : Integer := 0;

	begin

	   if matrix_rows /= matrix_cols then
	      raise Constraint_Error with "[SolveSystemHermSymm] source matrix must be square";
	   end if;

	   if rhs_rows /= matrix_rows then
	      raise Constraint_Error with "[SolveSystemHermSymm] matrix and rhs must have equal no. of rows";
	   end if;

	   SYSV ( UPLO  => 'L',
	          A     => matrix,
	          LDA   => matrix_rows,
	          N     => matrix_cols,
	          B     => rhs,
	          LDB   => rhs_rows,
	          NRHS  => rhs_cols,
	          IPIV  => pivots,
	          WORK  => short_vector,
	          LWORK => -1,
	          INFO  => return_code );

	   declare
	      work_vector_max : Constant Integer := Integer( short_vector(1).Re );
	      work_vector     : Complex_Vector (1 .. work_vector_max);
	   begin

	      SYSV ( UPLO  => 'L',
	             A     => matrix,
	             LDA   => matrix_rows,
	             N     => matrix_cols,
	             B     => rhs,
	             LDB   => rhs_rows,
	             NRHS  => rhs_cols,
	             IPIV  => pivots,
	             WORK  => work_vector,
	             LWORK => work_vector_max,
	             INFO  => return_code );

	   end;

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystemHermSymm] matrix may be singular";
	   end if;

	   return rhs;

	end SolveSystemHermSymm;

	procedure SolveSystem
	   (Solution   : out Real_Vector;
	    Source_mat : Real_Matrix;
	    Source_rhs : Real_Vector;
	    Size       : Integer)
	is

	   matrix        : Real_Matrix (1..size,1..size);
	   mixed_sol_rhs : Real_Matrix (1..size,1..1);
	   pivots        : Integer_Vector (1..size);

	   return_code   : Integer := 0;

	begin

	   if size > matrix'length(1) then
	      raise Constraint_Error with "[SolveSystem] no. rows in matrix less than size";
	   end if;

	   if size > matrix'length(2) then
	      raise Constraint_Error with "[SolveSystem] no. columns in matrix less than size";
	   end if;

	   for i in 1..size loop
	      mixed_sol_rhs (i,1) := Source_rhs (i);
	      for j in 1..size loop
	         matrix (i,j) := Source_mat (i,j);
	      end loop;
	   end loop;

	   GESV (A    => matrix,
	         LDA  => size,
	         N    => size,
	         IPIV => pivots,
	         B    => mixed_sol_rhs,
	         LDB  => size,
	         NRHS => 1,
	         INFO => return_code);

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystem] matrix may be singular";
	   end if;

	   for i in 1..size loop
	      solution (i) := mixed_sol_rhs (i,1);
	   end loop;

	end SolveSystem;

	procedure SolveSystem
	   (Solution   : out Complex_Vector;
	    Source_mat : Complex_Matrix;
	    Source_rhs : Complex_Vector;
	    Size       : Integer)
	is

	   matrix        : Complex_Matrix (1..size,1..size);
	   mixed_sol_rhs : Complex_Matrix (1..size,1..1);
	   pivots        : Integer_Vector (1..size);

	   return_code   : Integer := 0;

	begin

	   if size > matrix'length(1) then
	      raise Constraint_Error with "[SolveSystem] no. rows in matrix less than size";
	   end if;

	   if size > matrix'length(2) then
	      raise Constraint_Error with "[SolveSystem] no. columns in matrix less than size";
	   end if;

	   for i in 1..size loop
	      mixed_sol_rhs (i,1) := Source_rhs (i);
	      for j in 1..size loop
	         matrix (i,j) := Source_mat (i,j);
	      end loop;
	   end loop;

	   GESV (A    => matrix,
	         LDA  => size,
	         N    => size,
	         IPIV => pivots,
	         B    => mixed_sol_rhs,
	         LDB  => size,
	         NRHS => 1,
	         INFO => return_code);

	   if return_code /= 0 then
	      lapack_return_code := return_code;
	      raise Constraint_Error with "[SolveSystem] matrix may be singular";
	   end if;

	   for i in 1..size loop
	      solution (i) := mixed_sol_rhs (i,1);
	   end loop;

	end SolveSystem;

end Ada_Lapack.Extras;
