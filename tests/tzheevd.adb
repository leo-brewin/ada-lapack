with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tzheevd is
   
   type Real is digits 18;

   package Real_Arrays     is new Ada.Numerics.Generic_Real_Arrays (Real);
   package Complex_Types   is new Ada.Numerics.Generic_Complex_Types (Real);
   package Complex_Arrays  is new Ada.Numerics.Generic_Complex_Arrays (Real_Arrays, Complex_Types);
   
   package Real_Maths      is new Ada.Numerics.Generic_Elementary_Functions (Real);
   package Complex_Maths   is new Ada.Numerics.Generic_Complex_Elementary_Functions (Complex_Types);

   package Real_IO         is new Ada.Text_IO.Float_IO (Real);
   package Integer_IO      is new Ada.Text_IO.Integer_IO (Integer);
   package Complex_IO      is new Ada.Text_IO.Complex_IO (Complex_Types);

   package Lapack          is new Ada_Lapack(Real, Complex_Types, Real_Arrays, Complex_Arrays);
      
   use Lapack;
   
   use Real_Arrays;
   use Complex_Types;
   use Complex_Arrays;
   
   use Real_IO;
   use Integer_IO;
   use Complex_IO;
   
   use Real_Maths;
   use Complex_Maths;
   
   matrix : Complex_Matrix (1..4,1..4);
   
   matrix_rows : Integer := Matrix'Length (1);
   matrix_cols : Integer := Matrix'Length (2);
   
   eigenvalues      : Real_Vector (1..matrix_rows);
   eigenvalues_rows : Integer := matrix_rows;
   
   short_complex_vector : Complex_Vector (1..1);
   short_real_vector : Real_Vector (1..1);
   short_intg_vector : Integer_Vector (1..1);
   
   return_code  : Integer;

begin
   
   matrix:= (( (  3.40,  0.00), ( -2.36, -1.93), ( -4.68,  9.55), (  5.37, -1.23) ),
             ( ( -2.36,  1.93), (  6.94,  0.00), (  8.13, -1.47), (  2.07, -5.78) ),
             ( ( -4.68, -9.55), (  8.13,  1.47), ( -2.14,  0.00), (  4.68,  7.44) ),
             ( (  5.37,  1.23), (  2.07,  5.78), (  4.68, -7.44), ( -7.42,  0.00) ));

   HEEVD ( JOBZ   => 'V',
           UPLO   => 'L',
           A      => matrix,
           N      => matrix_cols,
           LDA    => matrix_rows,
           W      => eigenvalues,
           WORK   => short_complex_vector,
           LWORK  => -1,
			  RWORK  => short_real_vector,
			  LRWORK => -1,
           IWORK  => short_intg_vector,
           LIWORK => -1,
           INFO   => return_code );

   declare
      complex_work_max : Constant Integer := Integer( short_complex_vector(1).Re );
      real_work_max    : Constant Integer := Integer( short_real_vector(1) );
      intg_work_max    : Constant Integer := short_intg_vector(1);
      complex_work     : Complex_Vector (1 .. complex_work_max);
      real_work        : Real_Vector (1 .. real_work_max);
      intg_work        : Integer_Vector (1 .. intg_work_max);
   begin

      HEEVD ( JOBZ   => 'V',
              UPLO   => 'L',
              A      => matrix,
              N      => matrix_cols,
              LDA    => matrix_rows,
              W      => eigenvalues,
              WORK   => complex_work,
              LWORK  => complex_work_max,
				  RWORK  => real_work,
				  LRWORK => real_work_max,
              IWORK  => intg_work,
              LIWORK => intg_work_max,
              INFO   => return_code );

   end;

   if (return_code /= 0) then
      Put ("ZHEEVD failed, the return code was : ");
      Put ( return_code );
      New_line;
   else
      
      Put_line("The eigenvalues");
      for i in eigenvalues'range loop
         put(eigenvalues(i),3,4,0);
         put(" ");
      end loop;      
      new_line;
      
      new_line;
      
      Put_line("The eigenvectors");
      for i in 1..matrix_rows loop
         for j in 1..matrix_cols loop
            put(matrix(i,j),3,4,0);
            put(" ");
         end loop;
         new_line;
      end loop;
      
   end if;
         
end tzheevd;
