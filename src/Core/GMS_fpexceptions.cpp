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
#include <math.h>
#include <fenv.h>
#include "GMS_fpexceptions.h"
#include "GMS_config.h"
#include "GMS_indices.h"
#include "GMS_error_macros.h"





//
//	Implementation
//

//
//	Parametrised macros
//
#if !defined (FPEXCEPT_CHECK_COUNTER)
#define FPEXCEPT_CHECK_COUNTER(count) \
  do {	\
			if ((*count) > 0) (*count) = 0; \
  } while (0);
#endif

#if !defined (FPEXCEPT_CHECK_FENV_RETURN)
#define FPEXCEPT_CHECK_FENV_RETURN(err,msg)     \
	do{										    \
		if ((err) != 0) {                       \
		     PRINT_MESSAGE_VALUE((msg), (err))  \
			 ierr = -1;							\
			 return;							\
		 }										\
	} while (0);
#endif




void gms::math::is_denormalf32_present(const float * __restrict data,
				       const int32_t nx,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)


	if (process_slowly == true) {

		for (int32_t i = 0; i != nx; ++i) {
			if (fpclassify(data[i]) == FP_SUBNORMAL)
				*dnc += 1;
	    }
	}
	else {

		for (int32_t i = 0; i != nx; ++i) {
			if (fpclassify(data[i]) == FP_SUBNORMAL){
				*dnc = 1;
				return;
			}
		}
	}
}

void gms::math::is_denormalf32_present(const float * __restrict data,
				       const int32_t nx,
				       const int32_t ny,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)

	if (process_slowly == true) {
		for (int32_t i = 0; i != nx; ++i) {

			for (int32_t j = 0; j != ny; ++j) {
				if (fpclassify(data[I2D(i,j)]) == FP_SUBNORMAL)
					*dnc += 1;
			}
		}
	}
	else {
		for (int32_t i = 0; i != nx; ++i) {

			for (int32_t j = 0; j != ny; ++j) {
				if (fpclassify(data[I2D(i,j)]) == FP_SUBNORMAL){
					*dnc = 1;
					return;
				}
			}
		}
	}
}

void gms::math::is_denormalf32_present(const float * __restrict data,
				       const int32_t nx,
				       const int32_t ny,
				       const int32_t nz,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)

	if (process_slowly == true) {
		for (int32_t i = 0; i != nx; ++i) {
			for (int32_t j = 0; j != ny; ++j) {

				for (int32_t k = 0; k != nz; ++k) {
					if (fpclassify(data[I3D(i,j,k)]) == FP_SUBNORMAL)
						*dnc += 1;
				}
			}
		}
	}
	else {
		for (int32_t i = 0; i != nx; ++i) {
			for (int32_t j = 0; j != ny; ++j) {

				for (int_t k = 0; k != nz; ++k) {
					if (fpclassify(data[I3D(i,j,k)]) == FP_SUBNORMAL) {
						*dnc = 1;
						return;
					}
				}
			}
		}
	}
}

void gms::math::is_denormalf32_present(const float * __restrict data,
				       const int32_t nx,
				       const int32_t ny,
				       const int32_t nz,
				       const int32_t nw,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)

	if (process_slowly == true) {
		for (int32_t i = 0; i != nx; ++i) {
			for (int32_t j = 0; j != ny; ++j) {
				for (int32_t k = 0; k != nz; ++k) {

					for (int32_t l = 0; l != nw; ++l) {
						if (fpclassify(data[I4D(i,j,k,l)]) == FP_SUBNORMAL)
							*dnc += 1;
					}
				}
			}
		}
	}
	else {
		for (int32_t i = 0; i != nx; ++i) {
			for (int32_t j = 0; j != ny; ++j) {
				for (int32_t k = 0; k != nz; ++k) {

					for (int32_t l = 0; l != nw; ++l) {
						if (fpclassify(data[I4D(i,j,k,l)]) == FP_SUBNORMAL) {
							*dnc = 1;
							return;
						}
					}
				}
			}
		}
	}
}

void gms::math::is_denormalf64_present(const double * __restrict data,
				       const int32_t nx,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)

	if (process_slowly == true) {

		for (int32_t i = 0; i != nx; ++i) {
			if (fpclassify(data[i]) == FP_SUBNORMAL)
				*dnc += 1;
		}
	}
	else {

		for (int32_t i = 0; i != nx; ++i) {
			if (fpclassify(data[i]) == FP_SUBNORMAL) {
				*dnc = 1;
				return;
			}
		}
	}
}

void gms::math::is_denormalf64_present(const double * __restrict data,
				       const int32_t nx,
				       const int32_t ny,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)

	if (process_slowly == true) {
		for (int32_t i = 0; i != nx; ++i) {

			for (int32_t j = 0; j != ny; ++j) {
				if (fpclassify(data[I2D(i, j)]) == FP_SUBNORMAL)
					*dnc += 1;
			}
		}
	}
	else {
		for (int32_t i = 0; i != nx; ++i) {

			for (int32_t j = 0; j != ny; ++j) {
				if (fpclassify(data[I2D(i,j)]) == FP_SUBNORMAL) {
					*dnc = 1;
					return;
				}
			}
		}
	}
}

void gms::math::is_denormalf64_present(const double * __restrict data,
				       const int32_t nx,
				       const int32_t ny,
				       const int32_t nz,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)

	if (process_slowly == true) {
		for (int32_t i = 0; i != nx; ++i) {
			for (int32_t j = 0; j != ny; ++j) {

				for (int32_t k = 0; k != nz; ++k) {
					if (fpclassify(data[I3D(i, j, k)]) == FP_SUBNORMAL)
						*dnc += 1;
				}
			}
		}
	}
	else {
		for (int32_t i = 0; i != nx; ++i) {
			for (int32_t j = 0; j != ny; ++j) {

				for (int32_t k = 0; k != nz; ++k) {
					if (fpclassify(data[I3D(i,j,k)]) == FP_SUBNORMAL) {
						*dnc = 1;
						return;
					}
				}
			}
		}
	}
}

void gms::math::is_denormalf64_present(const double * __restrict data,
				       const int64_t nx,
				       const int64_t ny,
				       const int64_t nz,
				       const int64_t nw,
				       uint32_t * dnc,
				       const bool process_slowly) {
	FPEXCEPT_CHECK_COUNTER(dnc)
#define size ((nx)*(ny)*(nz)*(nw))
	if (process_slowly == true) {
		for (int64_t i = 0LL; i != nx; ++i) {
			for (int64_t j = 0LL; j != ny; ++j) {
				for (int64_t k = 0LL; k != nz; ++k) {
#if (size) < (L1_SIZE_DOUBLES)
#pragma prefetch data:0:8
#else
#pragma prefetch data:1:8
#endif
					for (int64_t l = 0LL; l != nw; ++l) {
						if (fpclassify(data[I4D(i, j, k, l)]) == FP_SUBNORMAL)
							*dnc += 1;
					}
				}
			}
		}
	}
	else {
		for (int64_t i = 0LL; i != nx; ++i) {
			for (int64_t j = 0LL; j != ny; ++j) {
				for (int64_t k = 0LL; k != nz; ++k) {
#if (size) < (L1_SIZE_DOUBLES)
#pragma prefetch data:0:8
#else
#pragma prefetch data:1:8
#endif
					for (int64_t l = 0LL; l != nw; ++l) {
						if (fpclassify(data[I4D(i, j, k, l)]) == FP_SUBNORMAL) {
							*dnc = 1;
							return;
						}
					}
				}
			}
		}
	}
}

void gms::math::is_abnormalf32(const float * __restrict data,
			       const int32_t nx,
			       uint32_t * count,
			       const bool process_slowly,
			       const uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {

		case 0: { // DENORMAL
					if (process_slowly == true) {
					
						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_SUBNORMAL)
								 *count += 1;
						}
					}
					else {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_SUBNORMAL) {
								*count = 1;
								return;
							}
						}
					}
		   }
		   break;

		case 1: {
					// NAN
					if (process_slowly == true) {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_NAN)
								*count += 1;
						}
					}
					else {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_NAN) {
								*count = 1;
								return;
							}
						}
					}
		    }
			break;

		case 2: {
					// INF
					if (process_slowly == true) {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_INFINITE)
								*count += 1;
						}
				 }
				 else {

					 for (int32_t i = 0; i != nx; ++i) {
						 if (fpclassify(data[i]) == FP_INFINITE) {
							 *count = 1;
							 return;
						 }
					 }
			    }
		  }
		  break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf32 -- domain 1D!!")
		}
	}
}

void gms::math::is_abnormalf32( const float * __restrict data,
				const int32_t nx,
				const int32_t ny,
				uint32_t * count,
				const bool process_slowly,
			        const uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {

		case 0: {
					// DENORMAL
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i,j)]) == FP_SUBNORMAL)
									*count += 1;
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_SUBNORMAL){
									*count = 1;
									return;
								}
							}
						}
					}
			}
			break;

		case 1: {
					// NAN
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i,j)]) == FP_NAN)
									*count += 1;
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_NAN) {
									*count = 1;
									return;
								}
							}
						}
					}
		    }
			break;

		case 2: {
					// INFINITE
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i,j)]) == FP_INFINITE)
									*count += 1;
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_INFINITE) {
									*count = 1;
									return;
								}
							}
						}
					}
		    }
			break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf32 -- domain 2D.")
		}
	}
}

void	gms::math::is_abnormalf32(const float * __restrict data,
				  const int32_t nx,
				  const int32_t ny,
				  const int32_t nz,
				  uint32_t * count,
				  const bool process_slowly,
				  const uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {

		case 0: {
					// DENORMAL
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i,j,k)]) == FP_SUBNORMAL)
										*count += 1;
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_SUBNORMAL) {
										*count = 1;
										return;
									}
								}
							}
						}
				  }
		    }
		    break;

		case 1: {
					 // NAN
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i,j,k)]) == FP_NAN)
										*count += 1;
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_NAN) {
										*count = 1;
										return;
									}
								}
							}
						}
					}
			 }
			 break;

		case 2: {
					// INFINITE
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i,j,k)]) == FP_INFINITE)
										*count += 1;
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_INFINITE) {
										*count = 1;
										return;
									}
								}
							}
						}
				}
		  }
		  break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf32 -- domain 3D.")
		}
	}
}

void gms::math::is_abnormalf32(const float * __restrict data,
			       const int32_t nx,
			       const int32_t ny,
			       const int32_t nz,
			       const int32_t nw,
			       uint32_t  * count,
			       const bool process_slowly,
			       const uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {
			
		case 0: {
					// SUBNORMAL
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i,j,k,l)]) == FP_SUBNORMAL)
											 *count += 1;
									}
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_SUBNORMAL) {
											*count = 1;
											return;
										}
									}
								}
							}
						}
					}
		     }
			 break;

		case 1: {
					// NAN
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i,j,k,l)]) == FP_NAN)
											*count += 1;
									}
								}
							}
						}
				    } 
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_NAN) {
											*count = 1;
											return;
										}
									}
								}
							}
						}
					}
			  }
			  break;

		case 2: {
					// INFINITE
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int_t k = 0; k != nz; ++k) {

									for (int_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i,j,k,l)]) == FP_INFINITE)
											*count += 1;
									}
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_INFINITE) {
											*count = 1;
											return;
										}
										
									}
								}
							}
						}
					}
		    }
			break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf32 -- domain 4D.")
		}
	}
}

void gms::math::is_abnormalf64(const double * __restrict data,
			       const int32_t nx,
			       uint32_t * count,
			       const bool process_slowly,
			       const uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {

		case 0: { // DENORMAL
					if (process_slowly == true) {
					
						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_SUBNORMAL)
								*count += 1;
						}
					}
					else {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_SUBNORMAL) {
								*count = 1;
								return;
							}
						}
					}
		}
		break;

		case 1: {
					// NAN
					if (process_slowly == true) {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_NAN)
								*count += 1;
						}
					}
					else {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_NAN) {
								*count = 1;
								return;
							}
						}
					}
		}
		break;

		case 2: {
					// INF
					if (process_slowly == true) {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_INFINITE)
								*count += 1;
						}
					}
					else {

						for (int32_t i = 0; i != nx; ++i) {
							if (fpclassify(data[i]) == FP_INFINITE) {
								*count = 1;
								return;
							}
						}
					}
		}
		break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf64 -- domain 1D!!")
		}
	}
}

void gms::math::is_abnormalf64(const double * __restrict data,
			       const int32_t nx,
			       const int32_t ny,
			       uint32_t * count,
			       const bool process_slowly,
			       uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {

		case 0: {
					// DENORMAL
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_SUBNORMAL)
									*count += 1;
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_SUBNORMAL){
									*count = 1;
									return;
								}
							}
						}
					}
		}
		break;

		case 1: {
					// NAN
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_NAN)
									*count += 1;
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_NAN) {
									*count = 1;
									return;
								}
							}
						}
					}
		}
		break;

		case 2: {
					// INFINITE
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_INFINITE)
									*count += 1;
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {

							for (int32_t j = 0; j != ny; ++j) {
								if (fpclassify(data[I2D(i, j)]) == FP_INFINITE) {
									*count = 1;
									return;
								}
							}
						}
					}
		}
		break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf64 -- domain 2D.")
		}
	}
}

void gms::math::is_abnormalf64(const double * __restrict data,
			       const int32_t nx,
			       const int32_t ny,
			       const int32_t nz,
			       uint32_t * count,
			       const bool process_slowly,
			       const uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {

		case 0: {
					// DENORMAL
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_SUBNORMAL)
										*count += 1;
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_SUBNORMAL) {
										*count = 1;
										return;
									}
								}
							}
						}
					}
		}
		break;

		case 1: {
					// NAN
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_NAN)
										*count += 1;
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_NAN) {
										*count = 1;
										return;
									}
								}
							}
						}
					}
		}
		break;

		case 2: {
					// INFINITE
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_INFINITE)
										*count += 1;
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {

								for (int32_t k = 0; k != nz; ++k) {
									if (fpclassify(data[I3D(i, j, k)]) == FP_INFINITE) {
										*count = 1;
										return;
									}
								}
							}
						}
					}
		}
		break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf64 -- domain 3D.")
		}
	}
}

void gms::math::is_abnormalf64(const double * __restrict data,
			       const int32_t nx,
			       const int32_t ny,
			       const int32_t nz,
			       const int32_t nw,
			       uint32_t * count,
			       const bool process_slowly,
			       const uint32_t option) {
	FPEXCEPT_CHECK_COUNTER(count)

		switch (option) {

		case 0: {
					// SUBNORMAL
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0LL; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_SUBNORMAL)
											*count += 1;
									}
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_SUBNORMAL) {
											*count = 1;
											return;
										}
									}
								}
							}
						}
					}
		}
		break;

		case 1: {
					// NAN
					if (process_slowly == true) {
						for (int32_t i = 0LL; i != nx; ++i) {
							for (int32_t j = 0LL; j != ny; ++j) {
								for (int32_t k = 0LL; k != nz; ++k) {

									for (int32_t l = 0LL; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_NAN)
											*count += 1;
									}
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0LL; j != ny; ++j) {
								for (int32_t k = 0LL; k != nz; ++k) {

									for (int32_t l = 0LL; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_NAN) {
											*count = 1;
											return;
										}
									}
								}
							}
						}
					}
		}
		break;

		case 2: {
					// INFINITE
					if (process_slowly == true) {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_INFINITE)
											*count += 1;
									}
								}
							}
						}
					}
					else {
						for (int32_t i = 0; i != nx; ++i) {
							for (int32_t j = 0; j != ny; ++j) {
								for (int32_t k = 0; k != nz; ++k) {

									for (int32_t l = 0; l != nw; ++l) {
										if (fpclassify(data[I4D(i, j, k, l)]) == FP_INFINITE) {
											*count = 1;
											return;
										}

									}
								}
							}
						}
					}
		}
		break;

		default: {
					 PRINT_ERROR_INFO("Invalid parameter to switch in: is_abnormalf64 -- domain 4D.")
		}
	}
}

int32_t gms::math::clear_fpexcepts(void) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_ALL_EXCEPT);
		if(err != 0) {
			PRINT_ERROR_VALUE("Function: 'feclearexcept' -- failed with an error: ",err)
				return (err);
		}
		return (err);
	}
#if (SILENCE_COMPILER) == 1
	return (err); // if somehow this branch is executed return err = -9999.
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::clear_fedenormal(void) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_DENORMAL);
		if(err != 0) {
			PRINT_ERROR_VALUE("Function: 'feclearexcpet(FE_DENORMAL)' -- failed with an error: ",err)
				return (err);
		}
		return (err);
	}
#if (SILENCE_COMPILER) == 1
	return (err); // if somehow this branch is executed return err = -9999.
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::clear_feinexact(void) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_INEXACT);
		if(err != 0) {
			PRINT_ERROR_VALUE("Function: 'feclearexcept(FE_INEXACT)' -- failed with an error: ",err)
				return (err);
		}
		return (err);
	}
#if (SILENCE_COMPILER) == 1
	return (err); // if somehow this branch is executed return err = -9999.
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::clear_feinvalid(void) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT){
		err = feclearexcept(FE_INVALID);
		if (err != 0) {
			PRINT_ERROR_VALUE("Function: 'feclearexcept(FE_INVALID)' -- failed with an error: ",err)
			return (err);
		}
		return (err);
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::clear_fedivbyzero(void) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT){
		err = feclearexcept(FE_DIVBYZERO);
		if (err != 0) {
			PRINT_ERROR_VALUE("Function: 'feclearexcept(FE_DIVBYZERO)' -- failed with an error: ",err)
			return (err);
		}
		return (err);
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::clear_feoverflow(void) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT){
		err = feclearexcept(FE_OVERFLOW);
		if(err != 0) {
			PRINT_ERROR_VALUE("Function: 'feclearexcept(FE_OVERFLOW)' -- failed with an error: ",err)
			return (err);
		}
		return (err);
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::clear_feunderflow(void) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = feclearexcept(FE_UNDERFLOW);
		if(err != 0) {
			PRINT_ERROR_VALUE("Function: 'feclearexcept(FE_UNDEFLOW)' -- failed with an error: ",err)
			return (err);
		}
		return (err);
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::test_feexcepts(const int32_t fex) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = fetestexcept(fex);
		if(FE_ALL_EXCEPT == err) {
			PRINT_MESSAGE("'fetestexcept' detected, that all fe-exceptions are set!!")
			return (err);
		}
		else {
			PRINT_MESSAGE("'fetestexcept' detected, that some fe-exceptions are set, but not all!!")
			return (err);
         }
	}
#if (SILENCE_COMPILER) == 1
	return (err); // if somehow this branch is executed return err = -9999.
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::test_feinvalid(const int32_t feinv) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = fetestexcept(feinv);
		if (FE_INVALID == err) {
			PRINT_MESSAGE_VALUE("'fetestexcept' detected  FE_INVALID exception, numeric value: ",err)
			return (err);
		}
		else {
			PRINT_MESSAGE_VALUE("'fetestexcept' failed to detected FE_INVALID exception, numeric value: ",err)
			return (err);
		}
	}
#if (SILENCE_COMPILER) == 1
	return (err); // if somehow this branch is executed return err = -9999.
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::test_feinexact(const int32_t feinex) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = fetestexcept(feinex);
		if (FE_INEXACT == err) {
			PRINT_MESSAGE_VALUE("'fetestexcept' detected FE_INEXACT exception, numeric value: ",err)
			return (err);
		}
		else {
			PRINT_MESSAGE_VALUE("'fetestexcept' failed to detected FE_INEXACT exception, numeric value: ",err)
			return (err);
		}
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::test_fedivbyzero(const int32_t fedbz) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = fetestexcept(fedbz);
		if(FE_DIVBYZERO == err) {
			PRINT_MESSAGE_VALUE("'fetestexcept' detected FE_DIVBYZERO exception, numeric value: ",err)
			return (err);
		}
		else {
			PRINT_MESSAGE_VALUE("'fetestxcept' failed to detect FE_DIVBYZERO exception, numeric value: ",err)
			return (err);
        }
	}
#if (SILENCE_COMPILER) == 1
	 return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::test_feunormal( const int32_t feden) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = fetestexcept(feden);
		if(FE_DENORMAL == err) {
			PRINT_MESSAGE_VALUE("fetestexcept' detected FE_DENORMAL exception, numeric value: ",err)
			return (err);
		}
		else {
			PRINT_MESSAGE_VALUE("fetestexcept' failed to detected FE_DENORMAL exception, numeric value: ",err)
			return (err);
		}
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::test_feoverflow(const int32_t feov) {
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = fetestexcept(feov);
		if(FE_OVERFLOW == err) {
			PRINT_MESSAGE_VALUE("fetestexcept' detected FE_OVERFLOW exception, numeric value: ",err)
			return (err);
		}
		else {
			PRINT_MESSAGE_VALUE("fetestexcept' failed to detect FE_OVERFLOW exception, numeric value: ",err)
			return (err);
		}
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

int32_t gms::math::test_feunderflow(const int32_t feun) {
	int32_t err{ -9999 };
#pragma STD FENV_ACCESS ON
#if defined (math_errhandling) && defined (MATH_ERREXCEPT)
	if (math_errhandling & MATH_ERREXCEPT) {
		err = fetestexcept(feun);
		if (FE_UNDERFLOW == err) {
			PRINT_MESSAGE_VALUE("fetestexcept' detected FE_UNDERFLOW exception, numeric value: ", err)
				return (err);
		}
		else {
			PRINT_MESSAGE_VALUE("fetestexcept' failed to detect FE_UNDERFLOW exception, numeric value: ", err)
				return (err);
		}
	}
#if (SILENCE_COMPILER) == 1
	return (err);
#endif
#else
#error "Undefined: 'math_errorhandling and MATH_ERREXCEPT' "
#endif
}

void gms::math::rise_fedenormal(const bool clear_prev, 
				const int32_t feden,
				int32_t & ierr) {
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(FE_DENORMAL == feden);
#else
	if(FE_DENORMAL != feden) {
		PRINT_ERROR_VALUE("rise_fedenormal: Invalid argument: ",feden)
		ierr = -1;
		return;
	}
#endif
	if (ierr < 0) ierr = 0;
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
	if (clear_prev == true) {
		if (test_feunormal(feden) == FE_DENORMAL) {
			err = clear_fedenormal();
			FPEXCEPT_CHECK_FENV_RETURN("'clear_fedenormal' -- failed with an value: ", err)
			
			err = feraiseexcept(feden);
			FPEXCEPT_CHECK_FENV_RETURN("feraiseexcept failed with an value: ", err)
			
			ierr = 0;
		}
		else {
			PRINT_MESSAGE(" 'test_feunormal -- failed!!'")
			ierr = -2;
			return;
		}
	}
	else {
		err = feraiseexcept(feden);
		FPEXCEPT_CHECK_FENV_RETURN("feraiseexcept failed with an value : ",err)
		
		ierr = 0;
	}	
}	

void GMS::math::rise_feinvalid(const bool clear_prev,
			       const int32_t feinv,
			       int32_t & ierr) {
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(FE_INVALID == feinv);
#else
	if(FE_INVALID != feinv) {
		PRINT_MESSAGE_VALUE(" 'rise_feinvalid: ' Invalid argument: ",feinv)
		ierr = -1;
		return;
	}
#endif
	if (ierr < 0) ierr = 0;
	int32_t err{-9999};
#pragma STD FENV_ACCESS ON
	if (clear_prev == true) {
		if (test_feinvalid(feinv) == FE_INVALID) {
			err = clear_feinvalid();
			FPEXCEPT_CHECK_FENV_RETURN(" 'clear_feinvalid' -- failed with a value: ", err)
			
			err = feraiseexcept(feinv);
			FPEXCEPT_CHECK_FENV_RETURN(" 'feraiseexcept'  -- failed with a value: ", err)
			ierr = 0;
		}
		else {
			PRINT_MESSAGE(" 'test_feinvalid ' -- failed")
			ierr = -2;
			return;
		}
	}
	else {
		err = feraiseexcept(feinv);
		FPEXCEPT_CHECK_FENV_RETURN(" 'feraiseexcept' --failed with a value : ",err)
		
		ierr = 0;
	}
}

void GMS::math::rise_feinexact(const bool clear_prev,
			       const int32_t feinex,
			       int32_t & ierr) {
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(FE_INEXACT == feinex);
#else
	if(FE_INEXACT != feinex) {
		PRINT_MESSAGE_VALUE(" 'rise_feinexact: ' Invalid argument: ",feinex)
		ierr = -1;
		return;
	}
#endif
	if (ierr < 0) ierr = 0;
	int32_t err = -9999;
#pragma STD FNV_ACCESS ON
	if (clear_prev == true) {
		if (test_feinexact(feinex) == FE_INEXACT){
			err = clear_feinexact();
			FPEXCEPT_CHECK_FENV_RETURN(" 'clear_feinexact: ' failed with a value : ", err)
			
			err = feraiseexcept(feinex);
			FPEXCEPT_CHECK_FENV_RETURN(" ' feraiseexcept: ' failed with a value:", err)
			ierr = 0;
		}
		else {
			PRINT_MESSAGE(" 'test_feinexact' -- failed!!")
			ierr = -2;
			return;
		}
	}
	else {
		err = feraiseexcept(feinex);
		FPEXCEPT_CHECK_FENV_RETURN(" 'feraiseexcept' --failed with a value : ", err)
		ierr = 0;
	}
}

void gms::math::rise_fedivbyzero(const bool clear_prev,
				 const int32_t fdbz,
				 int32_t & ierr) {
#if (GMS_DEBUG_ON) == 1
	_ASSERTE(FE_DIVBYZERO == fdbz);
#else
	if(FE_DIVBYZERO != fdbz) {
		PRINT_MESSAGE_VALUE(" 'rise_fedivbyzero: ' Invalid argument: ", err)
		ierr = -1;
		return;
	}
#endif
	if (ierr < 0) ierr = 0;
	int32_t err = -9999;
#pragma STD FENV_ACCESS ON
	if (clear_prev == true) {
		if (test_fedivbyzero(fdbz) == FE_DIVBYZERO) {
			err = clear_fedivbyzero();
			FPEXCEPT_CHECK_FENV_RETURN(" 'clear_fedivbyzero: ' failed with a value: ", err)

			err = feraiseexcept(fdbz);
			FPEXCEPT_CHECK_FENV_RETURN(" 'feraiseexcept: ' failed with a value: ", err)
			ierr = 0;
		}
		else {
			PRINT_MESSAGE( " 'test_fedivbyzero: ' -- failed!!")
			ierr = -2;
			return;
		}
	}
	else {
		err = feraiseexcept(fdbz);
		FPEXCEPT_CHECK_FENV_RETURN(" 'feraiseexcept: ' -- failed with a value: ", err)
		ierr = 0;
	}
}

 void gms::math::rise_feoverflow(const bool clear_prev,
				 const int32_t feov,
				 int32_t & ierr) {
#if (GMS_DEBUG_ON) == 1
	 _ASSERTE(FE_OVERFLOW == feov);
#else
	 if(FE_OVERFLOW != feov) {
		 PRINT_MESSAGE_VALUE(" 'rise_feoverflow: ' Invalid argument: ",feov)
		 ierr = -1;
		 return;
	 }
#endif
	 if (ierr < 0) ierr = 0;
	 int32_t err = -9999;
#pragma STD FNV_ACCESS ON
	 if (clear_prev == true) {
		 if (test_feoverflow(feov) == FE_OVERFLOW) {
			 err = clear_feoverflow();
			 FPEXCEPT_CHECK_FENV_RETURN(" 'clear_feoverlow' failed with a value: ",err)

			 err = feraiseexcept(feov);
			 FPEXCEPT_CHECK_FENV_RETURN(" 'feraiseexcept' failed with a value: ", err)
			 ierr = 0;
		 }
		 else {
			 PRINT_MESSAGE(" 'test_feoverflow' -- failed!!")
			 ierr = -2;
			 return;
		 }
	 }
	 else {
		 err = feraiseexcept(feov);
		 FPEXCEPT_CHECK_FENV_RETURN(" ' feraiseexcept ' -- failed with a value", err)
		 ierr = 0;
	 }
}

 void gms::math::rise_feundeflow(const bool clear_prev,
				 const int32_t feun,
				 int32_t & ierr) {
#if (GMS_DEBUG_ON) == 1
	 _ASSERTE(FE_UNDERFLOW == feun);
#else
	 if(FE_UNDERFLOW != feun){
		 PRINT_MESSAGE_VALUE(" ' rise_feunderflow: ' Invalid argument: ",feun)
		 ierr = -1;
		 return;
	 }
#endif
	 if (ierr < 0) ierr = 0;
	 int32_t err = -9999;
#pragma STD FENV_ACCESS ON
	 if (clear_prev == true) {
		 if (test_feunderflow(feun) == FE_UNDERFLOW) {
			 err = clear_feunderflow();
			 FPEXCEPT_CHECK_FENV_RETURN(" 'clear_feunderflow: ' failed with a value: ", err)

			 err = feraiseexcept(feun);
			 FPEXCEPT_CHECK_FENV_RETURN(" 'feraiseexcept: ' failed with a value: ", err)
			 ierr = 0;
		 }
		 else {
			 PRINT_MESSAGE(" 'test_feunderflow: ' -- failed!!")
			 ierr = -2;
			 return;
		 }
	 }
	 else {
		 err = feraiseexcept(feun);
		 FPEXCEPT_CHECK_FENV_RETURN(" ' feraiseexcept: ' failed with a value: ",err)
		 ierr = 0;
	 }
 }



