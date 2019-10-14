with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tdsysv is
   
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
   
   matrix : Real_Matrix (1..5,1..5);
   
   matrix_rows : Integer := Matrix'Length (1);
   matrix_cols : Integer := Matrix'Length (2);
   
   rhs : Real_Matrix (1..5,1..3);
   
   rhs_rows : Integer := rhs'Length (1);
   rhs_cols : Integer := rhs'Length (2);
   
   pivots   : Integer_Vector (1..matrix_rows);
   
   short_vector : Real_Vector (1..1);
   
   return_code  : Integer;

begin
   
   matrix:= (( -5.86,   3.99,  -5.93,  -2.82,   7.69 ),
             (  3.99,   4.46,   2.58,   4.42,   4.61 ),
             ( -5.93,   2.58,  -8.52,   8.57,   7.69 ),
             ( -2.82,   4.42,   8.57,   3.72,   8.07 ),
             (  7.69,   4.61,   7.69,   8.07,   9.83 ));

   rhs:= ((  1.32,  -6.33,  -8.77 ),
          (  2.22,   1.69,  -8.33 ),
          (  0.12,  -1.56,   9.54 ),
          ( -6.41,  -9.49,   9.56 ),
          (  6.33,  -3.67,   7.48 ));

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

   if (return_code /= 0) then

      Put ("DSYSV failed, the return code was : ");
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

end tdsysv;
