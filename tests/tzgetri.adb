with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tzgetri is
   
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
   
   matrix      : Complex_Matrix (1..4,1..4);
   matrix_copy : Complex_Matrix (1..4,1..4);
   
   matrix_rows : Integer := Matrix'Length (1);
   matrix_cols : Integer := Matrix'Length (2);
         
   pivots : Integer_Vector (1..matrix_rows);

   short_vector : Complex_Vector (1..1);
   
   return_code : Integer;

begin
   
   matrix :=( ((  1.23, -5.50), (  7.91, -5.38), ( -9.80, -4.86), ( -7.32,  7.57)),
              (( -2.14, -1.12), ( -9.92, -0.79), ( -9.18, -1.12), (  1.37,  0.43)),
              (( -4.30, -7.10), ( -6.47,  2.52), ( -6.51, -2.67), ( -5.86,  7.38)),
              ((  1.27,  7.29), (  8.90,  6.92), ( -8.82,  1.25), (  5.41,  5.37)) );
   
   matrix_copy := matrix;
   
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
               
   end;

   if (return_code /= 0) then

      Put_line ("The matrix is probably singular.");

   else
      
      Put_line ("Matrix");
      for i in matrix_copy'range(1) loop
         for j in matrix_copy'range(2) loop
            put(matrix_copy(i,j),3,3,0);
            put(" ");
         end loop;
         new_line;
      end loop;
      new_line;
         
      Put_line ("Inverse");
      for i in matrix'range(1) loop
         for j in matrix'range(2) loop
            put(matrix(i,j),3,3,0);
            put(" ");
         end loop;
         new_line;
      end loop;
         
   end if;
      
end tzgetri;
