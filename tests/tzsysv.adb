with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tzsysv is
   
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
   
   rhs : Complex_Matrix (1..4,1..2);
   
   rhs_rows : Integer := rhs'Length (1);
   rhs_cols : Integer := rhs'Length (2);
   
   pivots   : Integer_Vector (1..matrix_rows);
   
   short_vector : Complex_Vector (1..1);
   
   return_code  : Integer;

begin
   
   matrix:= (( (  9.99, -4.73), ( -5.68, -0.80), ( -8.94,  1.32), ( -9.42,  2.05) ),
             ( ( -5.68, -0.80), ( -8.01,  4.61), (  1.64, -6.29), (  6.79, -2.17) ),
             ( ( -8.94,  1.32), (  1.64, -6.29), (  9.04,  3.96), ( -4.51, -7.54) ),
             ( ( -9.42,  2.05), (  6.79, -2.17), ( -4.51, -7.54), (  0.40,  4.06) ));

   rhs:= ((  (  5.71, -1.20), (  2.84, -0.18) ),
          (  ( -7.70,  6.47), ( -8.29, -1.72) ),
          (  (  3.77, -7.40), ( -4.28, -8.25) ),
          (  ( -3.78,  0.33), ( -2.70, -0.39) ));

   SYSV ( UPLO  => 'U',
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

      SYSV ( UPLO  => 'U',
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

   if (return_code /= 0) then

      Put ("ZSYSV failed, the return code was : ");
      Put ( return_code );
      New_line;

   else
      
      Put_line ("Solution");
      for i in rhs'range(1) loop
         for j in rhs'range(2) loop
            put(" ");
            put(rhs(i,j),3,3,0);
         end loop;
         new_line;
      end loop;
         
      new_line;
      Put_line ("The matrix factorization");
      for i in matrix'range(1) loop
         for j in matrix'range(2) loop
            put(" ");
            put(matrix(i,j),3,3,0);
         end loop;
         new_line;
      end loop;
         
      new_line;
      Put_line ("The pivot indices");
      put("(");
      for i in pivots'range loop
         put(pivots(i),3);
      end loop;
      put_line(" )");
      
   end if;

end tzsysv;
