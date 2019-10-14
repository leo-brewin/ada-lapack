with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tdgesvd is
   
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
   
   function min(left,right : Integer) return Integer is
   begin
      if left < right
         then return left;
         else return right;
      end if;
   end min;
   
   one  : constant Real := 1.0e0;
   zero : constant Real := 0.0e0;
   
   matrix : Real_Matrix (1..6,1..5);
   
   matrix_rows : Integer := Matrix'Length (1);
   matrix_cols : Integer := Matrix'Length (2);
   
   num_singular : Integer;
   singular_values : Real_Vector (1..min(matrix_rows,matrix_cols));
   
   u_rows   : Integer := matrix_rows;
   u_cols   : Integer := matrix_rows;
   u_matrix : Real_Matrix (1..u_rows,1..u_cols);
   
   v_rows      : Integer := matrix_cols;
   v_cols      : Integer := matrix_cols;
   v_transpose : Real_Matrix (1..v_rows,1..v_cols);
   
   short_vector : Real_Vector (1..1);
   
   return_code  : Integer;
   
begin
   
   matrix :=( (  8.79,   9.93,   9.83,   5.45,   3.16),
              (  6.11,   6.91,   5.04,  -0.27,   7.98),
              ( -9.15,  -7.93,   4.86,   4.85,   3.01),
              (  9.57,   1.64,   8.83,   0.74,   5.80),
              ( -3.49,   4.02,   9.80,  10.00,   4.27),
              (  9.84,   0.15,  -8.99,  -6.02,  -5.31) );

   GESVD ( JOBU   => 'S',
           JOBVT  => 'S',
           M      => matrix_rows,
           N      => matrix_cols,
           A      => matrix,
           LDA    => matrix_rows,
           S      => singular_values,
           U      => u_matrix,
           LDU    => u_rows,
           VT     => v_transpose,
           LDVT   => v_rows,
           WORK   => short_vector,
           LWORK  => -1,
           INFO   => return_code );
         
   declare
      work_vector_rows : Integer := Integer( short_vector(1) );
      work_vector : Real_Vector (1 .. work_vector_rows);
   begin

      GESVD ( JOBU   => 'S',
              JOBVT  => 'S',
              M      => matrix_rows,
              N      => matrix_cols,
              A      => matrix,
              LDA    => matrix_rows,
              S      => singular_values,
              U      => u_matrix,
              LDU    => u_rows,
              VT     => v_transpose,
              LDVT   => v_rows,
              WORK   => work_vector,
              LWORK  => work_vector_rows,
              INFO   => return_code );
               
   end;
   
   if return_code > 0 then
      Put ("DGESVD failed, the return code was : ");
      Put ( return_code );
      New_line;
   else
   
      num_singular := min(matrix_rows,matrix_cols);
      
      put_line("Singular values");
      for i in 1..num_singular loop
         put(singular_values(i),3,4,0);
      end loop;
      new_line;
      new_line;
         
      put_line("Matrix U");
      for i in 1..matrix_rows loop
         for j in 1..num_singular loop
            put(u_matrix(i,j),3,4,0);
         end loop;
         new_line;
      end loop;
      new_line;
         
      put_line("Matrix V transpose");
      for i in 1..num_singular loop
         for j in 1..matrix_cols loop
            put(v_transpose(i,j),3,4,0);
         end loop;
         new_line;
      end loop;
      
   end if;

end tdgesvd;
