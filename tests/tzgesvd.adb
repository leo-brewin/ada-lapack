with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tzgesvd is
   
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
   
   matrix : Complex_Matrix (1..3,1..4);
   
   matrix_rows : Integer := Matrix'Length (1);
   matrix_cols : Integer := Matrix'Length (2);
   
   num_singular : Integer;
   singular_values : Real_Vector (1..min(matrix_rows,matrix_cols));
   
   u_rows   : Integer := matrix_rows;
   u_cols   : Integer := matrix_rows;
   u_matrix : Complex_Matrix (1..u_rows,1..u_cols);
   
   -- note: must declare v_hermitian as square even though bottom rows may not be used (when rows>cols)
      
   v_rows      : Integer := matrix_cols;
   v_cols      : Integer := matrix_cols;
   v_hermitian : Complex_Matrix (1..v_rows,1..v_cols);
   
   short_vector : Complex_Vector (1..1);
   real_work_vector : Real_Vector (1..5*min(matrix_rows,matrix_cols));
   
   return_code  : Integer;
   
begin
   
   matrix := ( ( ( 5.91, -5.69), ( 7.09,  2.72), ( 7.78, -4.06), ( -0.79, -7.21) ),
               ( (-3.15, -4.08), (-1.89,  3.27), ( 4.57, -2.07), ( -3.88, -3.30) ),
               ( (-4.89,  4.20), ( 4.10, -6.70), ( 3.28, -3.84), (  3.84,  1.19) ) );

   GESVD ( JOBU   => 'A',
           JOBVT  => 'A',
           M      => matrix_rows,
           N      => matrix_cols,
           A      => matrix,
           LDA    => matrix_rows,
           S      => singular_values,
           U      => u_matrix,
           LDU    => u_rows,
           VT     => v_hermitian,
           LDVT   => v_rows,
           WORK   => short_vector,
           LWORK  => -1,
           RWORK  => real_work_vector,
           INFO   => return_code );
        
   declare
      work_vector_rows : Integer := Integer( Re(short_vector(1)) );
      work_vector : Complex_Vector (1 .. work_vector_rows);
   begin
   
      GESVD ( JOBU   => 'A',
              JOBVT  => 'A',
              M      => matrix_rows,
              N      => matrix_cols,
              A      => matrix,
              LDA    => matrix_rows,
              S      => singular_values,
              U      => u_matrix,
              LDU    => u_rows,
              VT     => v_hermitian,
              LDVT   => v_rows,
              WORK   => work_vector,
              LWORK  => work_vector_rows,
              RWORK  => real_work_vector,
              INFO   => return_code );
               
   end;
   
   if return_code > 0 then
      Put ("ZGESVD failed, the return code was : ");
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
         
      put_line("Matrix V Hermitian transpose");
      for i in 1..num_singular loop
         for j in 1..matrix_cols loop
            put(v_hermitian(i,j),3,4,0);
         end loop;
         new_line;
      end loop;
         
      new_line;
      
   end if;

end tzgesvd;
