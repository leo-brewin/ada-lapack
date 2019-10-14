with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tdgeev is
   
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
   
   real_eigenvalues   : Real_Vector (1..matrix_rows);
   imag_eigenvalues   : Real_Vector (1..matrix_rows);
   
   eigenvectors_rows  : Integer := matrix_rows;
   eigenvectors_cols  : Integer := matrix_rows;

   left_eigenvectors  : Real_Matrix (1..eigenvectors_rows,1..eigenvectors_cols);
   right_eigenvectors : Real_Matrix (1..eigenvectors_rows,1..eigenvectors_cols);
   
   short_vector : Real_Vector (1..1);
   
   return_code  : Integer;

   i,j : integer;

begin
   
   matrix:= (( -1.01,   0.86,  -4.60,   3.31,  -4.81 ),
             (  3.98,   0.53,  -7.04,   5.29,   3.55 ),
             (  3.30,   8.26,  -3.89,   8.20,  -1.51 ),
             (  4.43,   4.96,  -7.66,  -7.33,   6.18 ),
             (  7.31,  -6.43,  -6.16,   2.47,   5.58 ));

   GEEV ( JOBVL   => 'V',
          JOBVR   => 'V',
          A       => matrix,
          LDA     => matrix_rows,
          N       => matrix_cols,
          WR      => real_eigenvalues,
          WI      => imag_eigenvalues,
          VL      => left_eigenvectors,
          VR      => right_eigenvectors,
          LDVL    => eigenvectors_rows,
          LDVR    => eigenvectors_rows,
          WORK    => short_vector,
          LWORK   => -1,
          INFO    => return_code );
   
   declare
      work_vector_rows : Integer := Integer( short_vector(1) );
      work_vector : Real_Vector (1 .. work_vector_rows);
   begin

      GEEV ( JOBVL   => 'V',
             JOBVR   => 'V',
             A       => matrix,
             N       => matrix_cols,
             LDA     => matrix_rows,
             WR      => real_eigenvalues,
             WI      => imag_eigenvalues,
             VL      => left_eigenvectors,
             VR      => right_eigenvectors,
             LDVL    => eigenvectors_rows,
             LDVR    => eigenvectors_rows,
             WORK    => work_vector,
             LWORK   => work_vector_rows,
             INFO    => return_code );

   end;
      
   if (return_code /= 0) then
      Put ("DGEEV failed, the return code was : ");
      Put ( return_code );
      New_line;
   else
      
      Put_line("The eigenvalues");
      for i in real_eigenvalues'range loop
         if imag_eigenvalues(i) /= 0.0e0 then
            put("(");
            put(real_eigenvalues(i),3,3,0);
            put(",");
            put(imag_eigenvalues(i),3,3,0);
            put(") ");
         else
            put(real_eigenvalues(i),3,3,0);
            put(" ");
         end if;
      end loop;      
      new_line;

      new_line;
      
      Put_line("The left eigenvectors");
      i := 1;
      while i <= eigenvectors_rows loop
         j := 1;
         while j <= eigenvectors_cols loop
            if imag_eigenvalues(j) /= 0.0e0 then
               put("(");
               put(left_eigenvectors(i,j),3,3,0);
               put(",");
               put(left_eigenvectors(i,j+1),3,3,0);
               put(") (");
               put(left_eigenvectors(i,j),3,3,0);
               put(",");
               put(-left_eigenvectors(i,j+1),3,3,0);
               put(") ");
               j := j+2;
            else
               put(left_eigenvectors(i,j),3,3,0);
               put(" ");
               j := j+1;
            end if;
         end loop;
         i := i+1;
         new_line;
      end loop;
         
      new_line;
      
      Put_line("The right eigenvectors");
      i := 1;
      while i <= eigenvectors_rows loop
         j := 1;
         while j <= eigenvectors_cols loop
            if imag_eigenvalues(j) /= 0.0e0 then
               put("(");
               put(right_eigenvectors(i,j),3,3,0);
               put(",");
               put(right_eigenvectors(i,j+1),3,3,0);
               put(") (");
               put(right_eigenvectors(i,j),3,3,0);
               put(",");
               put(-right_eigenvectors(i,j+1),3,3,0);
               put(") ");
               j := j+2;
            else
               put(right_eigenvectors(i,j),3,3,0);
               put(" ");
               j := j+1;
            end if;
         end loop;
         i := i+1;
         new_line;
      end loop;
      
   end if;
   
end tdgeev;
