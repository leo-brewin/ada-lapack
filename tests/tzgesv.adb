with Ada.Text_IO; 
with Ada.Text_IO.Complex_IO;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Complex_Arrays;
with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics.Generic_Complex_Elementary_Functions;
with Ada_Lapack;

use Ada.Text_IO;

procedure tzgesv is
   
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
   
   one  : constant Complex := (1.0e0,0.0e0);
   zero : constant Complex := (0.0e0,0.0e0);
   
   matrix      : Complex_Matrix (1..4,1..4);
   matrix_copy : Complex_Matrix (1..4,1..4);
   
   matrix_rows : Integer := Matrix'Length (1);
   matrix_cols : Integer := Matrix'Length (2);
   
   rhs : Complex_matrix (1..4,1..2);
   
   rhs_rows : Integer := Rhs'Length(1);
   rhs_cols : Integer := Rhs'Length(2);
   
   solution : Complex_Matrix := rhs;
   
   solution_rows : Integer := rhs_rows;
   solution_cols : Integer := rhs_cols;
   
   pivots : Integer_Vector (1..matrix_rows);
   
   return_code : Integer;

   Lower   : Complex_Matrix(1..matrix_rows,1..matrix_cols);
   Upper   : Complex_Matrix(1..matrix_rows,1..matrix_cols);
   Product : Complex_Matrix(1..matrix_rows,1..matrix_cols);
   
   swap    : integer;
   sum     : Complex;
   error   : Real;
   order   : Integer_Vector(1..matrix_rows);

begin
   
   matrix :=( ((  1.23, -5.50), (  7.91, -5.38), ( -9.80, -4.86), ( -7.32,  7.57)),
              (( -2.14, -1.12), ( -9.92, -0.79), ( -9.18, -1.12), (  1.37,  0.43)),
              (( -4.30, -7.10), ( -6.47,  2.52), ( -6.51, -2.67), ( -5.86,  7.38)),
              ((  1.27,  7.29), (  8.90,  6.92), ( -8.82,  1.25), (  5.41,  5.37)) );
  
   rhs := ( ((  8.33, -7.32), ( -6.11, -3.81)),
            (( -6.18, -4.80), (  0.14, -7.71)),
            (( -5.71, -2.80), (  1.41,  3.40)),
            (( -1.60,  3.08), (  8.54, -4.05)) );
          
   matrix_copy := matrix;
   
   GESV ( A       => matrix,
          LDA     => matrix_rows,
          N       => matrix_rows,
          IPIV    => pivots,
          B       => solution,
          LDB     => solution_rows,
          NRHS    => solution_cols,
          INFO    => return_code );
      
   if (return_code /= 0) then

      Put_line ("The matrix is probably singular.");

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
      Put_line ("The LU factorization");
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
      
      new_line;
      put_line("The original matrix");
      
      for i in 1..matrix_rows loop
         for j in 1..matrix_rows loop
            put(matrix_copy(i,j),3,3,0);
            put(" ");
         end loop;
         new_line;
      end loop;
         
      new_line;
      put_line("The reconstructed matrix");
      
      -- build L and U
      
      for i in 1..matrix_rows loop
         Lower(i,i) := One;
         for j in 1..i-1 loop
            Lower(i,j) := matrix(i,j);
            Lower(j,i) := Zero;
         end loop;
      end loop;
         
      for i in 1..matrix_rows loop
         Upper(i,i) := matrix(i,i);
         for j in i+1..matrix_rows loop
            Upper(i,j) := matrix(i,j);
            Upper(j,i) := Zero;
         end loop;
      end loop;
      
      -- construct an order vector from pivots
      
      for i in 1..matrix_rows loop
         order(i) := i;
      end loop;
         
      for i in 1..matrix_rows loop
         swap := order(i);
         order(i) := order(pivots(i));
         order(pivots(i)) := swap;
      end loop;
      
      -- compute P*L*U
      
      for i in 1..matrix_rows loop
         for j in 1..matrix_rows loop
            sum := Zero;
            for k in 1..matrix_rows loop
               sum := sum + Lower(i,k)*Upper(k,j);
            end loop;
            Product(order(i),j) := sum;
         end loop;
      end loop;
         
      -- the reconstructed matrix
      
      error := 0.0e0;
      for i in 1..matrix_rows loop
         for j in 1..matrix_rows loop
            put(Product(i,j),3,3,0);
            put(" ");
            error := error + modulus(Product(i,j) - matrix_copy(i,j));
         end loop;
         new_line;
      end loop;
         
      new_line;
      put("The residual error : ");
      put(error,2,4);
      new_line;
         
   end if;
      
end tzgesv;
