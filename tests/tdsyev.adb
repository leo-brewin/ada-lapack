with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tdsyev is
   
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
   
   eigenvalues      : Real_Vector (1..matrix_rows);
   eigenvalues_rows : Integer := matrix_rows;
   
   short_vector : Real_Vector (1..1);
   
   return_code  : Integer;

begin
   
   matrix:= ((  1.96,  -6.49,  -0.47,  -7.20,  -0.65 ),
             ( -6.49,   3.80,  -6.39,   1.50,  -6.34 ),
             ( -0.47,  -6.39,   4.17,  -1.51,   2.67 ),
             ( -7.20,   1.50,  -1.51,   5.70,   1.80 ),
             ( -0.65,  -6.34,   2.67,   1.80,  -7.10 ));

   SYEV ( JOBZ  => 'V',
          UPLO  => 'U',
          A     => matrix,
          N     => matrix_cols,
          LDA   => matrix_rows,
          W     => eigenvalues,
          WORK  => short_vector,
          LWORK => -1,
          INFO  => return_code );

   declare
      work_vector_rows : Integer := Integer( short_vector(1) );
      work_vector : Real_Vector (1 .. work_vector_rows);
   begin

      SYEV ( JOBZ  => 'V',
             UPLO  => 'U',
             A     => matrix,
             N     => matrix_cols,
             LDA   => matrix_rows,
             W     => eigenvalues,
             WORK  => work_vector,
             LWORK => work_vector_rows,
             INFO  => return_code );

   end;

   if (return_code /= 0) then
      Put ("DSYEV failed, the return code was : ");
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
         end loop;
         new_line;
      end loop;
      
   end if;

end tdsyev;
