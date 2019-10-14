with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tzheev is
   
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
   
   short_vector : Complex_Vector (1..1);
   real_work_vector : Real_Vector (1..3*matrix_rows-2);
   
   return_code  : Integer;

begin
   
   matrix:= (( (  9.14,  0.00), ( -4.37, -9.22), ( -1.98, -1.72), ( -8.96, -9.50) ),
             ( ( -4.37,  9.22), ( -3.35,  0.00), (  2.25, -9.51), (  2.57,  2.40) ),
             ( ( -1.98,  1.72), (  2.25,  9.51), ( -4.82,  0.00), ( -3.24,  2.04) ),
             ( ( -8.96,  9.50), (  2.57, -2.40), ( -3.24, -2.04), (  8.44,  0.00) ));

   HEEV ( JOBZ  => 'V',
          UPLO  => 'L',
          A     => matrix,
          N     => matrix_cols,
          LDA   => matrix_rows,
          W     => eigenvalues,
          WORK  => short_vector,
          LWORK => -1,
          RWORK => real_work_vector,
          INFO  => return_code );

   declare
      work_vector_rows : Integer := Integer( short_vector(1).Re );
      work_vector : Complex_Vector (1 .. work_vector_rows);
   begin

      HEEV ( JOBZ  => 'V',
             UPLO  => 'L',
             A     => matrix,
             N     => matrix_cols,
             LDA   => matrix_rows,
             W     => eigenvalues,
             WORK  => work_vector,
             LWORK => work_vector_rows,
             RWORK => real_work_vector,
             INFO  => return_code );

   end;

   if (return_code /= 0) then
      Put ("ZHEEV failed, the return code was : ");
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

end tzheev;
