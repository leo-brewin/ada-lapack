with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;

with Ada_Lapack;
with Ada_Lapack.Extras;

use Ada.Text_IO;

procedure example02 is
   
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
   package Extras          is new Lapack.Extras;
      
   use Lapack;
   use Extras;
   
   use Real_Arrays;
   use Complex_Types;
   use Complex_Arrays;
   
   use Real_IO;
   use Integer_IO;
   use Complex_IO;
   
   use Real_Maths;
   use Complex_Maths;
   
   matrix       : Complex_Matrix (1..4,1..4);
   eigenvectors : Complex_Matrix (1..4,1..4);
   eigenvalues  : Complex_Vector (1..4);

begin
   
   matrix:= ((( -3.84,  2.25), ( -8.94, -4.75), (  8.95, -6.53), ( -9.87,  4.82)),
             (( -0.66,  0.83), ( -4.40, -3.82), ( -3.50, -4.26), ( -3.15,  7.36)),
             (( -3.99, -4.73), ( -5.88, -6.60), ( -3.36, -0.40), ( -0.75,  5.23)),
             ((  7.74,  4.18), (  3.66, -7.53), (  2.58,  3.60), (  4.59,  5.41)));

   EigenSystem(matrix,eigenvalues,eigenvectors);
   
   Put_line("The eigenvalues");
   for i in eigenvalues'range loop
      put(eigenvalues(i),3,4,0);
      put(" ");
   end loop;      
   new_line;

   new_line;
         
   Put_line("The right eigenvectors");
   for i in eigenvectors'range(1) loop
      for j in eigenvectors'range(2) loop
         put(eigenvectors(i,j),3,4,0);
      put(" ");
      end loop;
      new_line;
   end loop;
         
end example02;
