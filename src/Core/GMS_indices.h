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
#ifndef __GMS_INDICES_H__
#define __GMS_INDICES_H__

/*
	Flat multidiensional array indexing macros
*/

// To be used with CUDA
/*
     When used in this for loop configuration -- traversing is column-order
     for(i ...
         for(j ....
     When used in this for loop configuration -- traversing is row-order
     for(j ...
        for(i ....
*/
#define INDEX2(isize,i,j) i + isize*j
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))


/*#define I2D(i,j,dim0)                              i+j*dim0
#define I3D(i,j,dim0,k,dim1)                       i+j*dim0+k*dim0*dim1
#define I4D(i,j,dim0,k,dim1,l,dim2)                i+j*dim0+k*dim0*dim1+l*dim0*dim1*dim2
#define I5D(i,j,dim0,k,dim1,l,dim2,m,dim3)         i+j*dim0+k*dim0*dim1+l*dim0*dim1*dim2+m*dim0*dim1*dim2*dim3
#define I6D(i,j,dim0,k,dim1,l,dim2,m,dim3,n,dim4)  i+j*dim0+k*dim0*dim1+l*dim0*dim1*dim2+m*dim0*dim1*dim2*dim3+n*dim0*dim1*dim2*dim3*dim4
*/

//More optimized version
/*
     When used in this for loop configuration -- traversing is column-order
     for(i ...
         for(j ....
     When used in this for loop configuration -- traversing is row-order
     for(j ...
        for(i ....
*/
#define Index2D(i,j,dim0)                                        i+j*dim0
#define Index3D(i,j,dim0,k,dim01)                                i+j*dim0+k*dim01
#define Index4D(i,j,dim0,k,dim01,l,dim012)                       i+j*dim0+k*dim01+l*dim012
#define Index5D(i,j,dim0,k,dim01,l,dim012,m,dim0123)             i+j*dim0+k*dim01+l*dim012+m*dim0123
#define Index6D(i,j,dim0,k,dim01,l,dim012,m,dim0123,n,dim01234)  i+j*dim0+k*dim01+l*dim012+m*dim0123+n*dim01234

/*
	@Warning:
				Macro parameters named idim,jdim,kdim,ldim
				must be present in the calling scope.
*/
#define Idx2D(i,j)       ((i) * jdim + (j))

#define Idx3D(i,j,k)     ((i) * jdim + ((j) * kdim + (k)))

#define Idx4D(i,j,k,l)   ((i) * idim + ((j) * jdim + ((k) * kdim + (l))))

#define Idx5D(i,j,k,l,m) ((i) * idim + ((j) * jdim + (k) * kdim + (l) * ldim + (m)))

/*
	@Warning:
				Macro parameters named nx,ny,nz must be
				present in the calling scope.
*/

#define I2D(i,j)    ((i) * ny + (j))

#define I3D(i,j,k)  ((i) * ny + ((j) * nz + (k)))

#define I4D(i,j,k,l)  ((i) * nx + ((j) * ny + ((k) * nz + ((l)))))

/*
	In this case dimension parameters 
	must be passed to parametrized macro.
				
*/

#define Ix2D(i,ny,j) i * ny + j

#define Ix3D(i,ny,j,nz,k) i * ny + j * nz + k

// Different indexing scheme.


/*
	@Warning:
				Macro parameters named idim,jdim,kdim,ldim
				must be present in the calling scope.
*/
#define IS2D(i,idim,j)   i * jdim + j

#define IS3D(i,idim,j,jdim,k) i * idim + (j * jdim + k)

#define IS4D(i,idim,j,jdim,k,kdim,l) i * idim + (j * jdim + (k * kdim + l))

#define IS5D(i,idim,j,jdim,k,kdim,l,ldim,m) i * idim + (j * jdim + (k * kdim + (l * ldim + m)))

// WRF indexing scheme
// m_jme and m_kme must be present in the calling scope.
#define Dim3(j,k,i) ((j) + (m_jme) * ((k) + (m_kme) * (i)))

/*
 /// 1-d indexing operation 

        /// @param[in] i index for dimension 0
      
        T& operator()(long i) {

            return data[i*_stride[0]];
        }

     

        /// 2-d indexing operation 

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
      
        T& operator()(long i, long j) {

            return data[i*_stride[0]+j*_stride[1]];
        }

      
        /// 3-d indexing operation 

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
       
        T& operator()(long i, long j, long k) {

            return data[i*_stride[0]+j*_stride[1]+k*_stride[2]];
        }

     

        /// 4-d indexing operation 

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
     
        T& operator()(long i, long j, long k, long l) {

            return data[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]];
        }

     

        /// 5-d indexing operation 

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @param[in] m index for dimension 4
      
        T& operator()(long i, long j, long k, long l, long m) {

            return data[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]+m*_stride[4]];
        }

      

        /// 6-d indexing operation 

        /// @param[in] i index for dimension 0
        /// @param[in] j index for dimension 1
        /// @param[in] k index for dimension 2
        /// @param[in] l index for dimension 3
        /// @param[in] m index for dimension 4
        /// @param[in] n index for dimension 5
        
        T& operator()(long i, long j, long k, long l, long m, long n) {

            return data[i*_stride[0]+j*_stride[1]+k*_stride[2]+
                           l*_stride[3]+m*_stride[4]+n*_stride[5]];
        }

       
*/



#endif /*__GMS_INDICES_H__*/
