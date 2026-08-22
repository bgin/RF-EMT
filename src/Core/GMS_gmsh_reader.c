/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "GMS_gmsh_reader.h"

/******************************************************************************/

char ch_cap ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_CAP capitalizes a single character.

  Discussion:

    This routine should be equivalent to the library "toupper" function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 1998

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the character to capitalize.

    Output, char CH_CAP, the capitalized character.
*/
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
/******************************************************************************/

int ch_eqi ( char ch1, char ch2 )

/******************************************************************************/
/*
  Purpose:

    CH_EQI is TRUE (1) if two characters are equal, disregarding case.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char CH1, CH2, the characters to compare.

    Output, int CH_EQI, is TRUE (1) if the two characters are equal,
    disregarding case and FALSE (0) otherwise.
*/
{
  int value;

  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     
  if ( ch1 == ch2 )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

int ch_to_digit ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_TO_DIGIT returns the integer value of a base 10 digit.

  Example:

     CH  DIGIT
    ---  -----
    '0'    0
    '1'    1
    ...  ...
    '9'    9
    ' '    0
    'X'   -1

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the decimal digit, '0' through '9' or blank are legal.

    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
    character was 'illegal', then DIGIT is -1.
*/
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
/******************************************************************************/

void gmsh_data_read ( char *gmsh_filename, int node_dim, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GMSH_DATA_READ reads data from a GMSH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, character *GMSH_FILENAME, the GMSH filename.

    Input, int NODE_DIM, the spatial dimension.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_X[NODE_DIM*NODE_NUM], the node coordinates.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
    the nodes that make up each element.
*/
{
  char *error;
  int i;
  int ierror;
  FILE *input;
  int j;
  int k;
  int length;
  int level;
  char text[255];
  char* text_pointer;
  double x;

  input = fopen ( gmsh_filename, "rt" );

  if ( ! input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GMSH_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open input file \"%s\"\n", gmsh_filename );
    exit ( 1 );
  }

  level = 0;
 
  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Nodes" ) )
      {
        level = 1;
        j = 0;
      }
    }
    else if ( level == 1 )
    {
      s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        s_to_i4 ( text_pointer, &length, &ierror );
        text_pointer = text_pointer + length;
        for ( i = 0; i < node_dim; i++ )
        {
          x = s_to_r8 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
          node_x[i+j*node_dim] = x;
        }
        j = j + 1;
      }
    }
  }
/*
  Now read element information.
*/
  level = 0;

  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Elements" ) )
      {
        level = 1;
        j = 0;
      }
    }
    else if ( level == 1 )
    {
      s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndElements" ) )
      {
        break;
      }
      else
      {
        for ( k = 1; k <= 5; k++ )
        {
          s_to_i4 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
        }
        for ( i = 0; i < element_order; i++ )
        {
          k = s_to_i4 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
          element_node[i+j*element_order] = k;
        }
        j = j + 1;
      }
    }
  }

  fclose ( input );

  return;
}
/******************************************************************************/

void gmsh_size_read ( char *gmsh_filename, int *node_num, int *node_dim,
  int *element_num, int *element_order )

/******************************************************************************/
/*
  Purpose:

    GMSH_SIZE_READ reads sizes from a GMSH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, character *GMSH_FILENAME, the GMSH filename.

    Output, int *NODE_NUM, the number of nodes.

    Output, int *NODE_DIM, the spatial dimension.

    Output, int *ELEMENT_NUM, the number of elements.

    Output, int *ELEMENT_ORDER, the order of the elements.
*/
{
  char *error;
  int ierror;
  FILE *input;
  int k;
  int length;
  int level;
  const double r8_big = 1.0E+30;
  char text[255];
  char* text_pointer;
  double x;
  double x_max;
  double x_min;
  double y;
  double y_max;
  double y_min;
  double z;
  double z_max;
  double z_min;

  *node_num = 0;
  *node_dim = 0;

  x_max = - r8_big;
  x_min = + r8_big;
  y_max = - r8_big;
  y_min = + r8_big;
  z_max = - r8_big;
  z_min = + r8_big;

  input = fopen ( gmsh_filename, "rt" );

  if ( ! input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GMSH_SIZE_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open input file \"%s\"\n", gmsh_filename );
    exit ( 1 );
  }

  level = 0;
 
  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Nodes" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      *node_num = s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        s_to_i4 ( text_pointer, &length, &ierror );
        text_pointer = text_pointer + length;
        x = s_to_r8 ( text_pointer, &length, &ierror );
        x_min = r8_min ( x_min, x );
        x_max = r8_max ( x_max, x );
        text_pointer = text_pointer + length;
        y = s_to_r8 ( text_pointer, &length, &ierror );
        y_min = r8_min ( y_min, y );
        y_max = r8_max ( y_max, y );
        text_pointer = text_pointer + length;
        z = s_to_r8 ( text_pointer, &length, &ierror);
        text_pointer = text_pointer + length;
        z_min = r8_min ( z_min, z );
        z_max = r8_max ( z_max, z );
      }
    }
  }
/*
  Make a very simple guess as to the dimensionality of the data.
*/
  *node_dim = 3;
  if ( z_max == z_min )
  {
    *node_dim = 2;
    if ( y_max == y_min )
    {
      *node_dim = 1;
    }
  }
/*
  Now read element information.
*/
  level = 0;

  for ( ; ; )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, input );

    if ( !error )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text_pointer, "$Elements" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      *element_num = s_to_i4 ( text_pointer, &length, &ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text_pointer, "$EndElements" ) )
      {
        break;
      }
      else
      {
        k = 0;
        for ( ; ; )
        {
          s_to_i4 ( text_pointer, &length, &ierror );
          text_pointer = text_pointer + length;
          if ( ierror != 0 )
          {
            break;
          }
          k = k + 1;
        }
        *element_order = k - 5;
        break;
      }
    }
  }

  fclose ( input );

  return;
}
/******************************************************************************/

int *gmsh_mesh2d_element_data_example ( int element_num, int element_order )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH2D_ELEMENT_DATA_EXAMPLE: element information for the example.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_ORDER, the order of the elements.

    Output, int GMSH_MESH2D_ELEMENT_DATA_EXAMPLE[ELEMENT_ORDER*ELEMENT_NUM], 
    the indices of the nodes that make up each element.
*/
{
  int *element_node;
  int element_node_save[3*24] = {
    1,  2,  6,
    7,  6,  2,
    2,  3,  7,
    8,  7,  3,
    3,  4,  8,
    9,  8,  4,
    4,  5,  9,
   10,  9,  5,
    6,  7, 11,
   12, 11,  7,
    7,  8, 12,
   13, 12,  8,
    8,  9, 13,
   14, 13,  9,
    9, 10, 14,
   15, 14, 10,
   11, 12, 16,
   17, 16, 12,
   12, 13, 17,
   18, 17, 13,
   16, 17, 19,
   20, 19, 17,
   17, 18, 20,
   21, 20, 18 };

  element_node = i4mat_copy_new ( element_order, element_num,
    element_node_save );

  return element_node;
}
/******************************************************************************/

void gmsh_mesh2d_element_size_example ( int *element_num, int *element_order )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH2D_ELEMENT_SIZE_EXAMPLE: element size information for the example.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

  Author:

    John Burkardt

  Parameters:

    Output, int *ELEMENT_NUM, the number of elements.

    Output, int *ELEMENT_ORDER, the order of the elements.
*/
{
  *element_num = 24;
  *element_order = 3;

  return;
}
/******************************************************************************/

double *gmsh_mesh2d_node_data_example ( int node_num, int node_dim )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH2D_NODE_DATA_EXAMPLE returns the node information for the example.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, int NODE_DIM, the spatial dimension.

    Output, double GMSH_MESH2D_NODE_DATA_EXAMPLE[NODE_DIM*NODE_NUM], 
    the nodal coordinates.
*/
{
  double *node_coord;
  double node_coord_save[2*21] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    3.0, 0.0,
    4.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 1.0,
    0.0, 2.0,
    1.0, 2.0,
    2.0, 2.0,
    3.0, 2.0,
    4.0, 2.0,
    0.0, 3.0,
    1.0, 3.0,
    2.0, 3.0,
    0.0, 4.0,
    1.0, 4.0,
    2.0, 4.0 };

  node_coord = r8mat_copy_new ( 2, 21, node_coord_save );

  return node_coord;
}
/******************************************************************************/

void gmsh_mesh2d_node_size_example ( int *node_num, int *node_dim )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH2D_NODE_SIZE_EXAMPLE: sizes of node information for the example.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

  Author:

    John Burkardt

  Parameters:

    Output, int *NODE_NUM, the number of nodes.

    Output, int *NODE_DIM, the spatial dimension.
*/
{
  *node_num = 21;
  *node_dim = 2;

  return;
}
/******************************************************************************/

void gmsh_mesh1d_write ( char *gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH1D_WRITE writes 1D mesh data as a Gmsh mesh file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2014

  Author:

    John Burkardt

  Reference:

    Christophe Geuzaine, Jean-Francois Remacle,
    Gmsh: a three-dimensional finite element mesh generator with
    built-in pre- and post-processing facilities,
    International Journal for Numerical Methods in Engineering,
    Volume 79, Number 11, pages 1309-1331, 2009.

  Parameters:

    Input, char *GMSH_FILENAME, the name of the Gmsh file.

    Input, int M, the spatial dimension.

    Input, inte NODE_NUM, the number of nodes.

    Input, double NODE_X[M*NODE_NUM], the node coordinates.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
    the nodes that make up each element.
*/
{
  int element;
  int element_type;
  FILE *gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
/*
  Force 1-based indexing.
*/
  mesh_base_one ( node_num, element_order, element_num, element_node );
/*
  Open the file.
*/
  gmsh = fopen ( gmsh_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( gmsh, "$MeshFormat\n" );
  fprintf ( gmsh, "2.2 0 8\n" );
  fprintf ( gmsh, "$EndMeshFormat\n" );

  fprintf ( gmsh, "$Nodes\n" );
  fprintf ( gmsh, "%d\n", node_num );
  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( gmsh, "  %d  %g  0.0  0.0\n", 
      node + 1, node_x[0+node*m] );
  }
  fprintf ( gmsh, "$EndNodes\n" );

  element_type = 1;

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, "$Elements\n" );
  fprintf ( gmsh, "%d\n", element_num );
  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( gmsh, "  %d  %d  %d  %d  %d",
      element + 1, element_type, tag_num, tag1, element + 1 );
    for ( i = 0; i < element_order; i++ )
    {
      fprintf ( gmsh, "  %d", element_node[i+element*element_order] );
    }
    fprintf ( gmsh, "\n" );
  }
  fprintf ( gmsh, "$EndElements\n" );

  fclose ( gmsh );

  return;
}
/******************************************************************************/

void gmsh_mesh2d_write ( char *gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH2D_WRITE writes 2D mesh data as a Gmsh mesh file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2014

  Author:

    John Burkardt

  Reference:

    Christophe Geuzaine, Jean-Francois Remacle,
    Gmsh: a three-dimensional finite element mesh generator with
    built-in pre- and post-processing facilities,
    International Journal for Numerical Methods in Engineering,
    Volume 79, Number 11, pages 1309-1331, 2009.

  Parameters:

    Input, char *GMSH_FILENAME, the name of the Gmsh file.

    Input, int M, the spatial dimension.

    Input, inte NODE_NUM, the number of nodes.

    Input, double NODE_X[M*NODE_NUM], the node coordinates.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
    the nodes that make up each element.
*/
{
  int element;
  int element_type;
  FILE *gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
/*
  Force 1-based indexing.
*/
  mesh_base_one ( node_num, element_order, element_num, element_node );
/*
  Open the file.
*/
  gmsh = fopen ( gmsh_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( gmsh, "$MeshFormat\n" );
  fprintf ( gmsh, "2.2 0 8\n" );
  fprintf ( gmsh, "$EndMeshFormat\n" );

  fprintf ( gmsh, "$Nodes\n" );
  fprintf ( gmsh, "%d\n", node_num );
  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( gmsh, "  %d  %g  %g  0.0\n", 
      node + 1, node_x[0+node*2], node_x[1+node*2] );
  }
  fprintf ( gmsh, "$EndNodes\n" );

  if ( element_order == 3 )
  {
    element_type = 2;
  }
  else if ( element_order == 6 )
  {
    element_type = 9;
  }

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, "$Elements\n" );
  fprintf ( gmsh, "%d\n", element_num );
  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( gmsh, "  %d  %d  %d  %d  %d",
      element + 1, element_type, tag_num, tag1, element + 1 );
    for ( i = 0; i < element_order; i++ )
    {
      fprintf ( gmsh, "  %d", element_node[i+element*element_order] );
    }
    fprintf ( gmsh, "\n" );
  }
  fprintf ( gmsh, "$EndElements\n" );

  fclose ( gmsh );

  return;
}
/******************************************************************************/

void gmsh_mesh3d_write ( char *gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH3D_WRITE writes 3D mesh data as a Gmsh mesh file.

  Discussion:

    The node ordering for the 20 node element is not standard.

    Assuming the vertices are A, B, C and D, Gmsh uses the following ordering:

    1:    a
    2:        b
    3:            c
    4:                d
    5: (2*a  +b        )/3
    6: (  a+2*b        )/3
    7: (    2*b+  c    )/3
    8: (      b+2*c    )/3
    9: (  a    +2*c    )/3
   10: (2*a    +  c    )/3
   11: (2*a        +  d)/3
   12: (  a        +2*d)/3
   13: (     b     +2*d)/3
   14: (   2*b     +  d)/3
   15: (       +  c+2*d)/3
   16: (       +2*c+  d)/3
   17: (  a+  b+  c    )/3
   18: (  a+  b    +  d)/3
   19: (      b+  c+  d)/3
   20: (  a+      c+  d)/3

    Leo Rebholz used the following ordering:

    1:    a
    2:        b
    3:            c
    4:                d
    5: (2*a  +b        )/3
    6: (2*a    +  c    )/3
    7: (  a+2*b        )/3
    8: (  a    +2*c    )/3
    9: (  a+  b+  c    )/3
   10: (    2*b+  c    )/3
   11: (      b+2*c    )/3
   12: (2*a        +  d)/3
   13: (   2*b     +  d)/3
   14: (       +2*c+  d)/3
   15: (  a+  b    +  d)/3
   16: (      b+  c+  d)/3
   17: (  a+      c+  d)/3
   18: (  a        +2*d)/3
   19: (     b     +2*d)/3
   20: (       +  c+2*d)/3

    Since the only 20 node data we have is from Leo, we will assume that
    all 20 node input data is in Leo's format, and needs to be converted
    to the Gmsh convention.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2014

  Author:

    John Burkardt

  Reference:

    Christophe Geuzaine, Jean-Francois Remacle,
    Gmsh: a three-dimensional finite element mesh generator with
    built-in pre- and post-processing facilities,
    International Journal for Numerical Methods in Engineering,
    Volume 79, Number 11, pages 1309-1331, 2009.

  Parameters:

    Input, char *GMSH_FILENAME, the name of the Gmsh file.

    Input, int M, the spatial dimension.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_X[M*NODE_NUM], the node coordinates.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
    the nodes that make up each element.
*/
{
  int element;
  int element_type;
  FILE *gmsh;
  int i;
  int i2;
  int leo_to_gmsh[20] = {
     0,  1,  2,  3,  4, 
     6,  9, 10,  7,  5, 
    11, 17, 18, 12, 19, 
    13,  8, 14, 15, 16 };
  int node;
  int tag_num;
  int tag1;
/*
  Force 1-based indexing.
*/
  mesh_base_one ( node_num, element_order, element_num, element_node );
/*
  Open the file.
*/
  gmsh = fopen ( gmsh_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( gmsh, "$MeshFormat\n" );
  fprintf ( gmsh, "2.2 0 8\n" );
  fprintf ( gmsh, "$EndMeshFormat\n" );

  fprintf ( gmsh, "$Nodes\n" );
  fprintf ( gmsh, "%d\n", node_num );
  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( gmsh, "%d  %g  %g  %g\n", 
      node + 1, node_x[0+node*3], node_x[1+node*3], node_x[2+node*3] );
  }
  fprintf ( gmsh, "$EndNodes\n" );

  if ( element_order == 4 )
  {
    element_type = 4;
  }
  else if ( element_order == 10 )
  {
    element_type = 11;
  }
  else if ( element_order == 20 )
  {
    element_type = 29;
  }

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, "$Elements\n" );
  fprintf ( gmsh, "%d\n", element_num );
  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( gmsh, "%d  %d  %d  %d  %d",
      element + 1, element_type, tag_num, tag1, element + 1 );
    if ( element_order == 20 )
    {
      for ( i = 0; i < element_order; i++ )
      {
        i2 = leo_to_gmsh[i];
        fprintf ( gmsh, "  %d", element_node[i2+element*element_order] );
      }
      fprintf ( gmsh, "\n" );
    }
    else
    {
      for ( i = 0; i < element_order; i++ )
      {
        fprintf ( gmsh, "  %d", element_node[i+element*element_order] );
      }
      fprintf ( gmsh, "\n" );
    }
  }
  fprintf ( gmsh, "$EndElements\n" );

  fclose ( gmsh );

  return;
}
/******************************************************************************/

int *i4mat_copy_new ( int m, int n, int a1[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_COPY_NEW copies an I4MAT to a "new" I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, int A1[M*N], the matrix to be copied.

    Output, int I4MAT_COPY_NEW[M*N], the copy of A1.
*/
{
  int *a2;
  int i;
  int j;

  a2 = ( int * ) malloc ( m * n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
/******************************************************************************/

void i4mat_transpose_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    fprintf ( stdout, "\n" );
/*
  For each row I in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Row: " );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "%6d  ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < jhi )
    {
      j2hi = n;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
/*
  Print out (up to INCX) entries in column J, that lie in the current strip.
*/
      fprintf ( stdout, "%5d: ", j - 1 );
      for ( i = i2lo; i <= i2hi; i++ )
      {
        fprintf ( stdout, "%6d  ", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void mesh_base_one ( int node_num, int element_order, int element_num,
  int element_node[] )

/******************************************************************************/
/*
  Purpose:

    MESH_BASE_ONE ensures that the element definition is one-based.

  Discussion:

    The ELEMENT_NODE array contains nodes indices that form elements.
    The convention for node indexing might start at 0 or at 1.

    If this function detects 0-based indexing, it converts to 1-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
    definitions.
*/
{
  int element;
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ONE:\n" );
    printf ( "  The element indexing appears to be 0-based!\n" );
    printf ( "  This will be converted to 1-based.\n" );
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] + 1;
      }
    }
  }
  else if ( node_min == 1 && node_max == node_num )
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ONE:\n" );
    printf ( "  The element indexing appears to be 1-based!\n" );
    printf ( "  No conversion is necessary.\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ONE - Warning!\n" );
    printf ( "  The element indexing is not of a recognized type.\n" );
    printf ( "  NODE_MIN = %d\n", node_min );
    printf ( "  NODE_MAX = %d\n", node_max );
    printf ( "  NODE_NUM = %d\n", node_num );
  }
  return;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Discussion:

    The C math library provides the function fmax() which should be preferred.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Discussion:

    The C math library provides the function fmin() which should be preferred.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
/******************************************************************************/

double *r8mat_copy_new ( int m, int n, double a1[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
*/
{
  double *a2;
  int i;
  int j;

  a2 = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return a2;
}
/******************************************************************************/

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "  %7d     ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "%5d:", j - 1 );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        fprintf ( stdout, "  %14g", a[(i-1)+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

int s_begin ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:

    S_BEGIN reports whether string 1 begins with string 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, char *S2, two strings.

    Output, int S_BEGIN, is true if S1 is the same as S2 up to
    the end of S2, and false otherwise.
*/
{
  int i;
  int n1;
  int n2;

  n1 = strlen ( s1 );
  n2 = strlen ( s2 );

  if ( n1 < n2 )
  {
    return 0;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return 0;
    }
  }
  return 1;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Discussion:

    It turns out that I also want to ignore the '\n' character!

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
    If S_LEN_TRIM is 0, then the string is entirely blank.
*/
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' && *t != '\n' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
/******************************************************************************/

int s_to_i4 ( char *s, int *last, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_I4 reads an I4 from a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a string to be examined.

    Output, int *LAST, the last character of S used to make IVAL.

    Output, int *ERROR is TRUE (1) if an error occurred and FALSE (0) otherwise.

    Output, int *S_TO_I4, the integer value read from the string.
    If the string is blank, then IVAL will be returned 0.
*/
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = 0;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  while ( *s ) 
  {
    c = s[i];
    i = i + 1;
/*
  Haven't read anything.
*/
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read the sign, expecting digits.
*/
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read at least one digit, expecting more.
*/
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
/*
  If we read all the characters in the string, see if we're OK.
*/
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = 1;
    *last = 0;
  }

  return ival;
}
/******************************************************************************/

double s_to_r8 ( char *s, int *lchar, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_R8 reads an R8 value from a string.

  Discussion:

    We have had some trouble with input of the form 1.0E-312.
    For now, let's assume anything less than 1.0E-20 is zero.

    This routine will read as many characters as possible until it reaches
    the end of the string, or encounters a character which cannot be
    part of the real number.

    Legal input is:

       1 blanks,
       2 '+' or '-' sign,
       2.5 spaces
       3 integer part,
       4 decimal point,
       5 fraction part,
       6 'E' or 'e' or 'D' or 'd', exponent marker,
       7 exponent sign,
       8 exponent integer part,
       9 exponent decimal point,
      10 exponent fraction part,
      11 blanks,
      12 final comma or semicolon.

    with most quantities optional.

  Example:

    S                 R

    '1'               1.0
    '     1   '       1.0
    '1A'              1.0
    '12,34,56'        12.0
    '  34 7'          34.0
    '-1E2ABCD'        -100.0
    '-1X2ABCD'        -1.0
    ' 2E-1'           0.2
    '23.45'           23.45
    '-4.2E+2'         -420.0
    '17d2'            1700.0
    '-14e-2'         -0.14
    'e2'              100.0
    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 June 2005

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string containing the
    data to be read.  Reading will begin at position 1 and
    terminate at the end of the string, or when no more
    characters can be read to form a legal real.  Blanks,
    commas, or other nonnumeric data will, in particular,
    cause the conversion to halt.

    Output, int *LCHAR, the number of characters read from
    the string to form the number, including any terminating
    characters such as a trailing comma or blanks.

    Output, int *ERROR, is TRUE (1) if an error occurred and FALSE (0)
    otherwise.

    Output, double S_TO_R8, the value that was read from the string.
*/
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = 0;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
/*
  Blank or TAB character.
*/
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
/*
  Comma.
*/
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
/*
  Minus sign.
*/
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Plus sign.
*/
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Decimal point.
*/
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Exponent marker.
*/
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Digit.
*/
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }
    }
/*
  Anything else is regarded as a terminator.
*/
    else
    {
      iterm = 1;
    }
/*
  If we haven't seen a terminator, and we haven't examined the
  entire string, go get the next character.
*/
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
/*
  If we haven't seen a terminator, and we have examined the
  entire string, then we're done, and LCHAR is equal to NCHAR.
*/
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
/*
  Number seems to have terminated.  Have we got a legal number?
  Not if we terminated in states 1, 2, 6 or 7!
*/
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = 1;
    return r;
  }
/*
  Number seems OK.  Form it.

  We have had some trouble with input of the form 1.0E-312.
  For now, let's assume anything less than 1.0E-20 is zero.
*/
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      if ( jsgn * jtop < -20 )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = pow ( ( double ) 10.0, ( double ) ( jsgn * jtop ) );
      }
    }
    else
    {
      if ( jsgn * jtop < -20 * jbot )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = jsgn * jtop;
        rexp = rexp / jbot;
        rexp = pow ( ( double ) 10.0, ( double ) rexp );
      }
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
