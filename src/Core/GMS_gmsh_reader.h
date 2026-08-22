
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
#ifndef __GMS_GMSH_READER_H__
#define __GMS_GMSH_READER_H__ 080120221349

#include <cstdint>
#include "GMS_config.h"


typedef struct __ATTR_ALIGN__(64) GMSH_Mesh {

       int32_t node_dim;            //NODE_DIM, the spatial dimension.
       int32_t node_num;            //NODE_NUM, the number of nodes.
       int32_t element_order;       //ELEMENT_ORDER
       int32_t element_num;         //ELEMENT_NUM, the number of elements
       double  * __restrict node_x;  //NODE_X[NODE_DIM*NODE_NUM], the node coordinates
       int32_t * __restrict element_node; //ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
                                          //the nodes that make up each element.
}GMSH_mesh;

char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void gmsh_data_read ( char *gmsh_filename, int node_dim, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] );
int *gmsh_mesh2d_element_data_example ( int element_num, int element_order );
void gmsh_mesh2d_element_size_example ( int * element_num, int *element_order );
double *gmsh_mesh2d_node_data_example ( int node_num, int node_dim );
void gmsh_mesh2d_node_size_example ( int *node_num, int *node_dim );
void gmsh_mesh1d_write ( char *gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh2d_write ( char *gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh3d_write ( char *gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_size_read ( char *gmsh_filename, int *node_num, int *node_dim,
  int *element_num, int *element_order );
int *i4mat_copy_new ( int m, int n, int a1[] );
void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r8mat_copy_new ( int m, int n, double a1[] );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
int s_begin ( char *s1, char *s2 );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
double s_to_r8 ( char *s, int *lchar, int *error );
void timestamp ( );


#endif /*__GMS_GMSH_READER_H__*/
